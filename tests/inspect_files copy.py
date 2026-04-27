import pickle

pkl_path = "/home/link/John_UMIs/bcwithqc/tests/pe_mini/pe_mini_threadsN_1/namepairidx_pe_mini_2_r1.pkl"

with open(pkl_path, "rb") as f:
    data = pickle.load(f)

print(data)