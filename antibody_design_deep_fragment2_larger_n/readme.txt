##evolutionary restrain  
## need first prepare the CDRs database, commands similar in prepare_all.bash，then:

bash  run_all_query_n.bash
##go to each replace_fragement_identify_point_?? folders, and change the "begin" and "ori_seq" according to given antibody residue id and CDRs
cd all_generate_single_m
bash cp_pkl.bash
python design_second_nn.py


##static potential

##static_eval and static_eval_cutoff10 are almost same, except that static_eval use cutoff 0.6nm as cutoff for antibody-antigen residue pair, while static_eval_cutoff10 use 1nm as cutoff
python caculate_interaction.py
bash run_score.bash

















