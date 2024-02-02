import pymysql
import pandas as pd
import base64
import os
import sys

def retrieve_ld(chr=None, rsid=None, population=None, rsid_col_name=None):
    print('当前文件名称: ',os.path.basename(__file__))
    print('当前函数名称: ',sys._getframe().f_code.co_name)
    try:
        host = base64.b64decode(
            'bG9jdXNjb21wYXJlLXVzLXdlc3QtMmEuY3ZvY3ViMzlubnJpLnVzLXdlc3QtMi5yZHMuYW1hem9uYXdzLmNvbQ==').decode()
        user = base64.b64decode('bG9jdXNjb21wYXJlcg==').decode()
        password = base64.b64decode('MTIzNDU2Nzg=').decode()
        database = base64.b64decode('bG9jdXNjb21wYXJl').decode()
        db = pymysql.connect(host=host, user=user, password=password, database=database)
        cursor = db.cursor()
        cursor.execute(f"select SNP_A, SNP_B, R2 from tkg_p3v5a_ld_chr{chr}_{population} where SNP_A = '{rsid}';")
        results_a = cursor.fetchall()

        cursor.execute(
            f"select SNP_B as SNP_A, SNP_A as SNP_B, R2 from tkg_p3v5a_ld_chr{chr}_{population} where SNP_B = '{rsid}';")
        results_b = cursor.fetchall()

        res = pd.concat([pd.DataFrame(list(results_a), columns=['target_snp', rsid_col_name, 'r2']),
                         pd.DataFrame(list(results_b), columns=['target_snp', rsid_col_name, 'r2']),
                         pd.DataFrame([{'target_snp': rsid, rsid_col_name: rsid, 'r2': 2}],
                                      columns=['target_snp', rsid_col_name, 'r2'])])

        cursor.execute(
            f"select EAS_AF, AMR_AF, AFR_AF, EUR_AF, SAS_AF, AF,ref, alt  from tkg_p3v5a_hg38 where rsid = '{rsid}';")
        results_rsid_af = cursor.fetchall()
        return res[[rsid_col_name, 'r2']], results_rsid_af[0]
    # once query r2 mysql db fail, only return target snp.
    except:
        return pd.DataFrame([{rsid_col_name: rsid, 'r2': 2}], columns=[rsid_col_name, 'r2']), (
        -1, -1, -1, -1, -1, -1, '', '')


# if __name__ == '__main__':
#     print(retrieve_ld(6, 'rs114907058', 'EUR', 'aa'))
