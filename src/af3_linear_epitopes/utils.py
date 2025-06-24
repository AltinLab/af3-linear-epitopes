import polars as pl
import hashlib
from mdaf3.FeatureExtraction import serial_apply, split_apply_combine
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from Bio import SeqIO
import polars as pl
import numpy as np


def fasta_to_polars(fasta_path: str, desc_as_name: bool = False) -> pl.DataFrame:
    """
    Read a FASTA file and convert it into a Polars DataFrame
    with columns ["name", "sequence"].

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file.

    Returns
    -------
    pl.DataFrame
        - "name": the sequence ID
        - "sequence": the full sequence string
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))

    if desc_as_name:
        names = [rec.description for rec in records]
    else:
        names = [rec.id for rec in records]
    seqs = [str(rec.seq) for rec in records]

    df = pl.DataFrame({"name": names, "seq": seqs})
    return df


def generate_job_name(df, cols, name="job_name"):
    df = df.with_columns(
        pl.concat_str(
            pl.concat_str(
                [
                    *[pl.col(colname) for colname in cols],
                ],
                ignore_nulls=True,
            )
            .map_elements(lambda x: hash_sequence(x, "md5"), return_dtype=pl.String)
            .alias(name),
        )
    )
    return df


def hash_sequence(seq: str, hash_type: str = "md5") -> str:
    """
    Hash a TCR sequence using the specified hash function.

    Args:
        tcr_seq (str): The TCR sequence string.
        hash_type (str): The hash function to use ('md5', 'sha1', 'sha256', etc.)

    Returns:
        str: The hexadecimal digest of the hashed sequence.
    """
    # Select the hash function
    if hash_type.lower() == "md5":
        h = hashlib.md5()
    elif hash_type.lower() == "sha1":
        h = hashlib.sha1()
    elif hash_type.lower() == "sha256":
        h = hashlib.sha256()
    else:
        raise ValueError("Unsupported hash type")

    # Encode the sequence and compute the hash
    h.update(seq.encode("utf-8"))
    return h.hexdigest()


def extract_sequence_from_ac(row):

    # uniprot is going to get mad at us for sending 60k requests
    # we don't care
    # just keep retrying until it works
    # https://stackoverflow.com/questions/23013220/max-retries-exceeded-with-url-in-requests
    session = requests.Session()
    retry = Retry(connect=10, backoff_factor=0.5)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("http://", adapter)
    session.mount("https://", adapter)

    ac = row["uniprotAC"]
    true_ac = ac

    params = {"fields": ["sequence"]}

    # base case: normal accession number in uniprot
    r = session.get(f"https://rest.uniprot.org/uniprotkb/{ac}", params=params)

    try:
        r.raise_for_status()
        j = r.json()

    # case 1: accession number not found. it's probably a cross reference identifier
    except Exception as e:
        print(f"Accession {ac} not valid, broadening search")

        params = {
            "query": ac,
            "fields": ["sequence"],
        }
        r = session.get(f"https://rest.uniprot.org/uniprotkb/search", params=params)
        r.raise_for_status()
        tmp_json = r.json()["results"]
        # make sure that it is a unique cross reference identifier
        if len(tmp_json) > 1:
            raise ValueError(f"Multiple query hits for accession {ac}")
        elif len(tmp_json) == 0:
            raise ValueError(f"No query hits for accession {ac}")
        else:
            j = tmp_json[0]

    # both cases have now returned the same json format
    # now, check if the entry was removed from unitProtKB and is now in uniparc
    if j["entryType"] == "Inactive":
        if j["inactiveReason"]["inactiveReasonType"] == "DEMERGED":
            # just grab first new ID
            true_ac = j["inactiveReason"]["mergeDemergeTo"][0]
            params = {"fields": ["sequence"]}

            r = session.get(
                f"https://rest.uniprot.org/uniprotkb/{true_ac}", params=params
            )
            r.raise_for_status()
            j = r.json()

        elif j["inactiveReason"]["inactiveReasonType"] == "DELETED":
            upi = j["extraAttributes"]["uniParcId"]

            params = {
                "upi": upi,
                "fields": ["sequence"],
            }

            r = session.get(
                "https://rest.uniprot.org/uniparc/%7Bupi%7D",
                params=params,
            )
            r.raise_for_status()
            j = r.json()
            true_ac = upi

        else:
            raise ValueError(
                f"Unkown reason {j["inactiveReason"]["inactiveReasonType"]} for {ac}"
            )

    try:
        row["sequence"] = j["sequence"]["value"]
    except KeyError:
        raise ValueError(f"No sequence found for accession {ac}")

    # since some accessions are removed
    # replace it with something that's currently working
    row["liveAC"] = true_ac

    return row


def set_30mer_indices_to_true(row_struct):
    """
    Sets specified indices in a boolean mask to True.

    Args:
        row_struct: A dictionary-like object representing a row, with keys
                    "boolmask" and "indices".

    Returns:
        A Polars Series containing the modified boolean mask.
    """
    boolmask = np.array(row_struct["boolmask"])
    indices = row_struct["indices"]
    for idx in indices:
        boolmask[idx : idx + 30] = True
    return boolmask.tolist()
