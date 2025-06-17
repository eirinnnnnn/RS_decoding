import subprocess
import os
from pathlib import Path
from tqdm import tqdm

def run_test(t1, t0, num=1000, build_cmd="./rs_decoder"):
    # Step 1: Generate input and answer files
    subprocess.run(["python3", "gentest.py", "-n", str(num), "-e", str(t1), "-E", str(t0)], check=True)

    # Step 2: Run decoder binary
    subprocess.run([build_cmd], check=True)

    # Step 3: Compare output.txt with answer.txt
    with open("output.txt") as f_out, open("answer.txt") as f_ans:
        out_lines = f_out.readlines()
        ans_lines = f_ans.readlines()
    
    assert len(out_lines) == len(ans_lines)
    # correct = 0
    # for out, ans in zip(out_lines, ans_lines):
    #     if out.strip() == ans.strip():
    #         correct += 1
    # return correct / len(ans_lines)
    false_accepts = 0
    for out, ans in zip(out_lines, ans_lines):
        if "*" not in out:
            false_accepts += 1  # should have failed, but decoded something
    return 1 - (false_accepts / len(ans_lines))  # detection rate


def main():
    Path("results").mkdir(exist_ok=True)
    # Assume rs_decoder is already compiled
    build_cmd = "./rs_decoder"

    print("Experiment I: fix t0=21 + 2·t1 = 21")
    with open("results/test1.csv", "w") as f:
        f.write("t0,t1,success_rate\n")
        # for a in range (0, 22):
        for t0 in tqdm(range(22, 32)):
            t1 = 0
            rate = run_test(t1=t1, t0=t0, build_cmd=build_cmd)
            f.write(f"{t0},{t1},{rate:.3f}\n")

    # print("Experiment II: fix t0 = 0, vary t1")
    # with open("results/test2.csv", "w") as f:
    #     f.write("t0,t1,success_rate\n")
    #     for t1 in tqdm(range(13, 20)):
    #         rate = run_test(t1=t1, t0=0, build_cmd=build_cmd)
    #         f.write(f"0,{t1},{rate:.3f}\n")

    # print("Experiment III: fix t1 = 0, vary t0")
    # with open("results/test3.csv", "w") as f:
    #     f.write("t0,t1,success_rate\n")
    #     for t0 in tqdm(range(25, 41)):
    #         rate = run_test(t1=0, t0=t0, build_cmd=build_cmd)
    #         f.write(f"{t0},0,{rate:.3f}\n")
    
    # print("Experiment IV: fix t0=19,20,21, t1=3")
    # with open("results/test4.csv", "w") as f:
    #     f.write("t0,t1,success_rate\n")
    #     for t0 in tqdm(range(16, 22)):
    #         t1 = 3 
    #         if t0 < 0: continue
    #         rate = run_test(t1=t1, t0=t0, build_cmd=build_cmd)
    #         f.write(f"{t0},{t1},{rate:.3f}\n")

if __name__ == "__main__":
    main()
