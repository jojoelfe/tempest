import pandas as pd
import sqlite3

def get_tm_results_from_database(database):
    con = sqlite3.connect(database)
    # I want:
    #   - ResID
    #   - JobID
    #   - Image
    #   - Template
    #   - Num matches
    #   - median score
    #   - max score
    tm_results = pd.read_sql_query(f"""
    SELECT 
        TEMPLATE_MATCH_ID, 
      
        IMAGE_ASSETS.NAME AS IMAGE_NAME,
        VOLUME_ASSETS.NAME AS REFERENCE_NAME

    FROM TEMPLATE_MATCH_LIST
    
    INNER JOIN IMAGE_ASSETS ON IMAGE_ASSETS.IMAGE_ASSET_ID = TEMPLATE_MATCH_LIST.IMAGE_ASSET_ID
    INNER JOIN VOLUME_ASSETS ON VOLUME_ASSETS.VOLUME_ASSET_ID = TEMPLATE_MATCH_LIST.REFERENCE_VOLUME_ASSET_ID""",con)
    
    tm_results['NUM_MATCHES'], tm_results['AVG_SCORE'], tm_results['MAX_SCORE']  = zip(*tm_results.apply(lambda row :
    pd.read_sql_query(f"SELECT COUNT(*), AVG(PEAK_HEIGHT), MAX(PEAK_HEIGHT) FROM TEMPLATE_MATCH_PEAK_LIST_{row['TEMPLATE_MATCH_ID']}",con).values[0],axis=1))
    con.close()
    return(tm_results)


if __name__ == "__main__":
    print(get_tm_results_from_database("C:\\Users\\jojot\\template_matching_chimerax\\patch_motion.db"))
