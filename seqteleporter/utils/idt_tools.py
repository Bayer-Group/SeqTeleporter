from __future__ import print_function
from base64 import b64encode
import json
from urllib import request, parse
import requests
from typing import Dict, List, Union, Type

from seqteleporter.utils.idt_config import IdtCredentials


def validate_seq_complexity_idt(dna: List[Dict[str, str]], idt_credentials: Type[IdtCredentials],
                                product_type: str) -> bool:
    print(f'Running IDT {product_type} complexity screener...')
    # divide into 100 sequences per batch (this is IDT API limitation)
    validated = True
    for batch_num in range(0, len(dna) // 100 + 1):
        dna_batch = dna[batch_num * 100:(batch_num + 1) * 100]
        complexity_result = idt_complexity_screener(dna_batch, idt_credentials, product_type)
        if type(complexity_result) == str:
            print(complexity_result)
            return False
        else:
            for idx, res in enumerate(complexity_result):
                if len(res) > 0:
                    validated = False
                    print(f'Detected problem with {dna[idx]["Name"]}.')
                    print(f'Sequence: {dna[idx]["Sequence"]}.')
                    for problem in res:
                        print(problem['DisplayText'])
                        # print(problem['ActualValue'])
                        # print(problem['ForwardLocations'])
    return validated


def idt_complexity_screener(dna: List[Dict[str, str]], idt_credentials: Type[IdtCredentials],
                            product_type: str) -> Union[List[List[dict]], str]:
    """
    See: https://eu.idtdna.com/restapi/swagger/ui/index#/Complexity
    :param dna: example input_draft: [{"Name": "stuffer","Sequence": "ATGCGTAGATT"}]
    :param product_type: "Gblock", "GblockHifi", "EBlock", "Gene"
    :return:
    """
    client_id, client_secret = idt_credentials.CLIENT_ID, idt_credentials.CLIENT_SECRET
    idt_username, idt_password = idt_credentials.IDT_USERNAME, idt_credentials.IDT_PASSWORD

    token = get_access_token(client_id, client_secret, idt_username, idt_password)
    url = ''.join(['https://eu.idtdna.com/restapi/v1/Complexities/Screen', product_type, 'Sequences'])
    max_seq_per_request = 100
    result_list: list = []
    if len(dna) > max_seq_per_request:
        for i in range(0, len(dna) // max_seq_per_request):
            dna_chunk = dna[max_seq_per_request * i:max_seq_per_request * (i + 1)]
            response = requests.post(url,
                                     headers={'Authorization': " ".join(["Bearer", token]),
                                              'Content-Type': 'application/json',
                                              'Accept': 'application/json'},
                                     data=str(dna_chunk))

            if response.status_code != 200:
                # print(response.text)
                return response.text
            result_list = result_list + json.loads(response.text)
    return result_list


def get_access_token(client_id, client_secret, idt_username, idt_password):
    """
    Create the HTTP request, transmit it, and then parse the response for the
    access token.

    The body_dict will also contain the fields "expires_in" that provides the
    time window the token is valid for (in seconds) and "token_type".
    """

    # Construct the HTTP request
    authorization_string = b64encode(bytes(client_id + ":" + client_secret, "utf-8")).decode()
    request_headers = {"Content-Type": "application/x-www-form-urlencoded",
                       "Authorization": "Basic " + authorization_string}

    data_dict = {"grant_type": "password",
                 "scope": "test",
                 "username": idt_username,
                 "password": idt_password}
    request_data = parse.urlencode(data_dict).encode()

    post_request = request.Request("https://eu.idtdna.com/Identityserver/connect/token",
                                   data=request_data,
                                   headers=request_headers,
                                   method="POST")

    # Transmit the HTTP request and get HTTP response
    response = request.urlopen(post_request)

    # Process the HTTP response for the desired data
    body = response.read().decode()

    # Error and return the response from the endpoint if there was a problem
    if response.status != 200:
        raise RuntimeError("Request failed with error code:" + response.status + "\nBody:\n" + body)

    body_dict = json.loads(body)
    return body_dict["access_token"]


if __name__ == "__main__":

    seq1 = 'AGGTAACCCAGTTCTTAGCACACATCCGTTTTCTCTATGACCACGCTCGATGTCGATCGCCTCAATTAGCGGACTTGTGTGCGTTAATGTGCTCCGTTGGGT' \
           'TGCCCCCAAGAAGTCGCCAAGAATTCATCGTAAGTGACCTGCGCTGTGGGCGACTATAGTGGTTTCTTCTGACGCAGCCTACCTGCGGTTTGAATTGCTGGT' \
           'CACATCGCTCTTCGACTTGCGGATGAAAACGCTTGCAGGCGAGCTTTACACTCCTAATATAAGCGGGCGAACCTGGCACGAGATACACTCTCGTTACCCAGT' \
           'GTGGTTGTGCTGCTGTACAAGACCATACGCAACTTAGGATGCGGGGTTTTTCCAATGCCCTACAACGGGGTCGCTATATGGACTTTCCAAGGCGAGCCCCGT' \
           'TGATTTTATAAGTAGCGGCCTAAGTGAATGCGACGTAGCGGCGTAAGGAGGAACTATTTGTCTACGCACGCTCACCGTATACTCCAGATCCG'
    seq2 = 'GCACACATCCGTTTTCTCTATGACCACGCTCGATGCACACATCCGTTTTCTCTATGACCACGCGTCGATCGCCTCAATTAGCGGACTTGTGAAAAAAAAAAAA'

    # dna = [{"Name": "stuffer", "Sequence": seq}]
    dna = [{"Name": f"seq{idx}", "Sequence": seq}
           for idx, seq in enumerate([seq1, seq2])]
    product_type = 'Gblock'
    idt_credentials = IdtCredentials
    complexity_result = idt_complexity_screener(dna, idt_credentials, product_type)
    if type(complexity_result) == str:
        print(complexity_result)
    else:
        for res in complexity_result:
            if len(res) > 0:
                for problem in res:
                    print(problem['DisplayText'])
                    print(problem['ActualValue'])
                    print(problem['ForwardLocations'])
