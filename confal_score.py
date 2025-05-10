import sys
import copy
import pandas as pd



def calculate_confal_score(df_row, mean_used_to_assign_df_row, cluster): 
    """
    observed_angles: list of 7 torsion angles
    cluster_center: list of mean torsion angles from conformer class
    std_devs: list of standard deviations of conformer class angles
    """

    angles_used_to_assign = ['d1', 'e1', 'z1', 'a2', 'b2', 'g2', 'd2']
    angles_stds_used_to_assign = ['d_1_std', 'e_1_std', 'z_1_std', 'a_2_std', 'b_2_std', 'g_2_std', 'd_1_std']
    
    mean_df_cluster = mean_used_to_assign_df_row[mean_used_to_assign_df_row['cluster'] == cluster]

    observed_angles = df_row[angles_used_to_assign].astype(float).tolist()
    cluster_center = mean_df_cluster[angles_used_to_assign].astype(float).iloc[0].tolist()
    std_devs = mean_df_cluster[angles_stds_used_to_assign].astype(float).iloc[0].tolist()

    score = 0

    for obs, center, std in zip(observed_angles, cluster_center, std_devs):
        delta = min(abs(obs - center), 360 - abs(obs - center))  # handle circular angles
        z = delta / std
        score += z**2
        
    return [cluster, max(0, 100 - score)]  # crude approximation





def assign_confal_score(df_row, mean_used_to_assign_df_row):

    confal_score_with_cluster = []
    for each_cluster in range(len(mean_used_to_assign_df_row)):
        confal_score_with_cluster.append(calculate_confal_score(df_row, mean_used_to_assign_df_row, each_cluster))
    
    sorted_confal_score_with_cluster = sorted(confal_score_with_cluster, key=lambda x: x[-1], reverse=True)

    return sorted_confal_score_with_cluster[0]






if __name__ == '__main__':

    rna_mean_df = pd.read_csv('/Users/savan/my_work/working/pro_rna/assignment/conf_assign/free_rna_conformers_present_mean.csv')
    rnp_mean_df = pd.read_csv('/Users/savan/my_work/working/pro_rna/assignment/conf_assign/rnp_conformers_present_mean.csv')
    
    if sys.argv[2] == 'rna':
        mean_df = rna_mean_df
    elif sys.argv[2] == 'rnp':
        mean_df = rnp_mean_df
    else:
        KeyError('Can assign either free RNA or RNA in RNP')


    to_assign_df = pd.read_csv(sys.argv[1], compression="gzip", header=0)

    to_assign_df[['cluster', 'confal_score']] = to_assign_df.apply(
    lambda row: assign_confal_score(row, mean_df), 
    axis=1, 
    result_type='expand'
)

    final_confal_score_df = copy.deepcopy(to_assign_df)


    final_confal_score_df.to_csv('assigned_cluster_with_confal_score.csv', index=False)

