import pymysql
import pandas as pd
import base64


def retrieve_ld(chr=None, rsid=None, population=None, rsid_col_name=None):
    try:
        host = base64.b64decode(
            'HOST').decode()
        user = base64.b64decode('USER').decode()
        password = base64.b64decode('PWD').decode()
        database = base64.b64decode('DB').decode()
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
        return res[[rsid_col_name, 'r2']]
    # once query r2 mysql db fail, only return target snp.
    except:
        return pd.DataFrame([{rsid_col_name: rsid, 'r2': 2}], columns=[rsid_col_name, 'r2'])


# if __name__ == '__main__':
#     print(retrieve_ld(16.0,'rs8055867','EUR','aa'))