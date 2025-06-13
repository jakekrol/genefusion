#!/usr/bin/env python3
import requests
import argparse

parser = argparse.ArgumentParser(
    description="Download sample table for a given COSMIC fusion ID."
)
parser.add_argument(
    "-i",
    "--fusion_id",
    type=str,
    required=True,
    help="COSMIC fusion ID to download sample table for."
)
args = parser.parse_args()


def download_fusion_sample_table(fusion_id, output="fusion_samples.tsv"):
    url = (
        f"https://cancer.sanger.ac.uk/cosmic/fusion/samples"
        f"?id={fusion_id}&iDisplayLength=10000&export=tsv"
    )

    headers = {
        "User-Agent": "Mozilla/5.0",
    }

    # go to devtools > network 
    # click the samples download button once
    # look for the request with sample table endpoint
    # example file
    # samples?id=1271&export=json&sEcho=1&iColumns=5&sColumns=&iDisplayStart=0&iDisplayLength=30&mDataProp_0=0&sSearch_0=&bRegex_0=false&bSearchable_0=true&bSortable_0=true&mDataProp_1=1&sSearch_1=&bRegex_1=false&bSearchable_1=true&bSortable_1=true&mDataProp_2=2&sSearch_2=&bRegex_2=false&bSearchable_2=true&bSortable_2=true&mDataProp_3=3&sSearch_3=&bRegex_3=false&bSearchable_3=true&bSortable_3=true&mDataProp_4=4&sSearch_4=&bRegex_4=false&bSearchable_4=true&bSortable_4=true&sSearch=&bRegex=false&iSortCol_0=0&sSortDir_0=asc&iSortingCols=1
    cookies = {
        "__Host-next-auth.csrf-token": "33c8db07913af7adf3433783538924f4dd7df8767be9a22793c6a078f6dd001d|82f1df2c4140d447eb1785ad4956b5960781ec18a92cbcf655bb70e602201ce2",
        "__Secure-next-auth.callback-url": "https://cancer.sanger.ac.uk/cosmic/download/cosmic",
        "__Secure-next-auth.session-token": "eyJhbGciOiJka... (truncated for clarity) ...WmFNn7Q.7_iPWpRKnj6OpERZpi2l7w",
        "CookiePolicy": "seen-e1-t1",
        "cosmic_session": "67958073039094375005298887527057060",
        "Pagesmith": "{\"z\":\"n\",\"a\":\"e\"}",
        "redirect_to": "/fusion/summary?id=1271"
    }

    response = requests.get(url, headers=headers, cookies=cookies)

    if response.ok:
        with open(output, "w", encoding="utf-8") as f:
            f.write(response.text)
        print(f"✅ Downloaded sample table for fusion ID {fusion_id} to: {output}")
    else:
        print(f"❌ Request failed: {response.status_code}")
        print(response.text)

# Example use
args.fusion_id = args.fusion_id.removeprefix("COSF")
download_fusion_sample_table(f"{args.fusion_id}", output=f"fusion_samples_{args.fusion_id}.tsv")
