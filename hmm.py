import numpy as np
import copy
import sys
import os, os.path

B_pseudocount = 0.04
A_pseudocount = 0.016

def filtering(old_seqs, old_seq_len, theta):

    filtered_seqs = copy.deepcopy(old_seqs)

    i_prim = 0
    for i in range(old_seq_len):       # moving on columns
        sum = 0
        for j in range(n):              # moving on seqs
            if old_seqs[j][i] == '-':
                sum += 1
        if sum > theta:
            filtered_seqs = np.delete(filtered_seqs, i_prim, 1)
        else:
            i_prim += 1
    return filtered_seqs

def DNA_or_Prot(seqs, n, T):

    for i in range(n):
        for j in range(T):
            if seqs[i][j] != 'A' and seqs[i][j] != 'C' and seqs[i][j] != 'G' \
                                 and seqs[i][j] != 'T' and seqs[i][j] != '-':
                return False
    return True

def convert_charToNum(seq_type, c):
    if seq_type:
        if c == 'A':
            return 0
        elif c == 'T':
            return 1
        elif c == 'G':
            return 2
        elif c == 'C':
            return 3
        elif c == '-':
            return 4
    else:
        if c == 'A':
            return 0
        elif c == 'R':
            return 1
        elif c == 'N':
            return 2
        elif c == 'D':
            return 3
        elif c == 'C':
            return 4
        elif c == 'Q':
            return 5
        elif c == 'E':
            return 6
        elif c == 'G':
            return 7
        elif c == 'H':
            return 8
        elif c == 'I':
            return 9
        elif c == 'L':
            return 10
        elif c == 'K':
            return 11
        elif c == 'M':
            return 12
        elif c == 'F':
            return 13
        elif c == 'P':
            return 14
        elif c == 'S':
            return 15
        elif c == 'T':
            return 16
        elif c == 'W':
            return 17
        elif c == 'Y':
            return 18
        elif c == 'V':
            return 19
        elif c == '-':
            return 20

def convert_NumToChar(seq_type, num):
    if seq_type:
        if num == 0:
            return 'A'
        elif num == 1:
            return 'T'
        elif num == 2:
            return 'G'
        elif num == 3:
            return 'C'
        elif num == 4:
            return '-'
    else:
        if num == 0:
            return 'A'
        elif num == 1:
            return 'R'
        elif num == 2:
            return 'N'
        elif num == 3:
            return 'D'
        elif num == 4:
            return 'C'
        elif num == 5:
            return 'Q'
        elif num == 6:
            return 'E'
        elif num == 7:
            return 'G'
        elif num == 8:
            return 'H'
        elif num == 9:
            return 'I'
        elif num == 10:
            return 'L'
        elif num == 11:
            return 'K'
        elif num == 12:
            return 'M'
        elif num == 13:
            return 'F'
        elif num == 14:
            return 'P'
        elif num == 15:
            return 'S'
        elif num == 16:
            return 'T'
        elif num == 17:
            return 'W'
        elif num == 18:
            return 'Y'
        elif num == 19:
            return 'V'
        elif num == 20:
            return '-'

def legals(T, num_of_observations):

    legal_pi = np.zeros(3*T+3, dtype=float)
    legal_A = np.zeros(shape=(3*T+3, 3*T+3), dtype=float)
    legal_B = np.zeros(shape=(3*T+3, num_of_observations), dtype=float)

    legal_pi[0] = 1
    for i in range(T):
        for j in range(num_of_observations-1):
            legal_B[3*i+1][j] = 1  # I
            legal_B[3*i+2][j] = 1   # M
        legal_B[3 * i + 3][num_of_observations - 1] = 1  # D

    for j in range(num_of_observations - 1):
        legal_B[3*T+1][j] = 1   # last I

    legal_B[0][num_of_observations-1] = 1  # S
    legal_B[3*T+2][num_of_observations-1] = 1  # E
    legal_A[0][1] = 1  # S->I0
    legal_A[0][2] = 1  # S->M1
    legal_A[0][3] = 1  # S->D1

    for i in range(T):

        legal_A[3*i+1][3*i+1] = 1  # I0->I0 ...
        legal_A[3*i+1][3*i+2] = 1  # I0->M1 ...

        legal_A[3*i+2][3*i+4] = 1  # M1->I1 ...
        legal_A[3*i+2][3*i+5] = 1  # M1->M2 ...
        if i != T-1:  # There's no transition from last M to a D state!
            legal_A[3*i+2][3*i+6] = 1  # M1->D2 ...

        legal_A[3*i+3][3*i+5] = 1  # D1->M2 ...
        if i != T-1:  # There's no transition from last D to another D state!
            legal_A[3*i+3][3*i+6] = 1  # D1->D2 ...

    legal_A[3*T+1][3*T+1] = 1  # last I->last I
    legal_A[3*T+1][3*T+2] = 1  # last I->E

    return legal_A, legal_B, legal_pi

def initialization(seqs, A, B, pi, T, n, num_of_observations):

    for i in range(0, T):

        # B matrix for Match/Mismatch states
        for j in range(n):
            B[3*i+2][int(seqs[j][i])] += 1

        zero_num = 0
        for j in range(num_of_observations-1):
            if B[3*i+2][j] == 0:
                zero_num += 1

        for j in range(num_of_observations-1):
            if B[3*i+2][j] == 0:
                B[3*i+2][j] = B_pseudocount
            else:
                B[3*i+2][j] = B[3*i+2][j] * (1 - zero_num * B_pseudocount) / n

        # B matrix for Insertion states
        for j in range(num_of_observations-1):
            B[3*i+1][j] = 1 / (num_of_observations-1)

        # B matrix for Deletion states
        B[3*i+3][num_of_observations-1] = 1

    for j in range(num_of_observations-1):        # last I state
        B[3*T+1][j] = 1 / (num_of_observations-1)

    B[0][num_of_observations - 1] = 1               # Start state
    B[3*T+2][num_of_observations - 1] = 1       # End state
    A[0][1] = 1/3                    # S->I0
    A[0][2] = 1/3                    # S->M1
    A[0][3] = 1/3                    # S->D1

    for i in range(T):

        # A matrix for Insertion states
        A[3*i+1][3*i+1] = 0.77           # I0->I0 ...
        A[3*i+1][3*i+2] = 0.23           # I0->M1 ...

        # A matrix for Match/Mismatch states
        A[3*i+2][3*i+4] = 0.025         # M1->I1 ...
        A[3*i+2][3*i+5] = 0.95          # M1->M2 ...
        if i != T-1:                    # There's no transition from last M to a D state!
            A[3*i+2][3*i+6] = 0.025     # M1->D2 ...
        else:
            A[3*i+2][3*i+4] = 0.5       # last M->last I
            A[3*i+2][3*i+5] = 0.5       # last M->E

        # A matrix for Deletion states
        A[3*i+3][3*i+5] = 0.33          # D1->M2 ...
        if i != T-1:                    # There's no transition from last D to another D state!
            A[3*i+3][3*i+6] = 0.67      # D1->D2 ...
        else:
            A[3*i+3][3*i+5] = 1         # last D->E

    A[3*T+1][3*T+1] = 0.77              # last I->last I
    A[3*T+1][3*T+2] = 0.23              # last I->E

    pi[0] = 1

    return A, B, pi

def create_alpha(seq, A, B, pi, T):

    alpha = np.zeros(shape=(3*T+3, T), dtype=float)
    for i in range(0, 3*T+3):
        alpha[i][0] = pi[i]*B[i][int(seq[0])]
    for t in range(1, T):
        for j in range(0, 3*T+3):
            sum = 0
            for i in range(3*T+3):
                sum += alpha[i][t-1] * A[i][j]
            alpha[j][t] = sum * B[j][int(seq[t])]

    return alpha

def create_beta(seq, A, B, T):

    beta = np.zeros(shape=(3*T+3, T), dtype=float)

    for i in range(0, 3*T+3):
        beta[i][T-1] = 1

    for t in range(T-2, -1, -1):
        for i in range(0, 3*T+3):
            for j in range(3*T+3):
                beta[i][t] += A[i][j] * B[j][int(seq[t+1])] * beta[j][t+1]

    return beta

def create_sai(seq, alpha, beta, A, B, T):     # ξ

    sai = np.zeros(shape=(3*T+3, 3*T+3, T), dtype=float)

    for t in range(0, T-1):
        d_sum = 0
        for i in range(3*T+3):
            for j in range(3*T+3):
                sai[i][j][t] = alpha[i][t] * A[i][j] * B[j][int(seq[t+1])] * beta[j][t+1]
                d_sum += sai[i][j][t]
        for i in range(3*T+3):
            for j in range(3*T+3):
                sai[i][j][t] /= d_sum

    return sai

def Baum_Welch(seqs, A, B, pi, T, n, num_of_observations):

    alpha = np.zeros(shape=(n, 3*T+3, T), dtype=float)
    beta = np.zeros(shape=(n, 3*T+3, T), dtype=float)
    sai = np.zeros(shape=(n, 3*T+3, 3*T+3, T), dtype=float)

    legal_A, legal_B, legal_pi = legals(T, num_of_observations)

    for iteration in range(20):

        for s in range(n):
            alpha[s] = create_alpha(seqs[s], A, B, pi, T)
            beta[s] = create_beta(seqs[s], A, B, T)
            sai[s] = create_sai(seqs[s], alpha[s], beta[s], A, B, T)

        for i in range(3*T+3):
            # updating A
            for j in range(3*T+3):
                if legal_A[i][j] == 1:
                    temp = A[i][j]
                    A[i][j] = 0
                    d_sum = 0
                    for r in range(n):
                        for t in range(0, T - 1):
                            A[i][j] += sai[r][i][j][t]
                            for k in range(3 * T + 3):
                                d_sum += sai[r][i][k][t]
                    if d_sum == 0:
                        A[i][j] = temp
                        continue
                    A[i][j] /= d_sum

            # updating B
            for k in range(num_of_observations):
                if legal_B[i][k] == 1:
                    temp = B[i][k]
                    B[i][k] = 0
                    d_sum = 0
                    for r in range(n):
                        for t in range(0, T):
                            sum = 0
                            for z in range(3 * T + 3):
                                sum += sai[r][i][z][t]
                            d_sum += sum
                            if int(seqs[r][t]) == k:
                                B[i][k] += sum
                    if d_sum == 0:
                        B[i][k] = temp
                        continue
                    B[i][k] /= d_sum

        for i in range(3*T+3):
            if i != 3*T+2:
                A_sum = 0
                for j in range(3*T+3):
                    A[i][j] = max(A_pseudocount, A[i][j])
                    A[i][j] *= legal_A[i][j]
                    A_sum += A[i][j]
                for j in range(3*T+3):
                    A[i][j] /= A_sum

            B_sum = 0
            for j in range(num_of_observations):
                B[i][j] = max(B_pseudocount, B[i][j])
                B[i][j] *= legal_B[i][j]
                B_sum += B[i][j]
            for j in range(num_of_observations):
                B[i][j] /= B_sum

    return A, B, pi

def get_delta(l):
    return float(l[3])

def Viterbi2(seq, A, B, pi, T, T2, seq_type):

    viterbi_list1 = []

    delta = np.zeros(shape=(3*T+3, T2), dtype=float)  # δ

    for i in range(0, 3*T+3):
        delta[i][0] = pi[i] * B[i][int(seq[0])]

    for t in range(1, T2):
        for j in range(0, 3*T+3):
            maxx = 0
            for i in range(3*T+3):
                if delta[i][t - 1] * A[i][j] > maxx:
                    maxx = delta[i][t-1] * A[i][j]
            if j%3 == 0 and j != 0:
                delta[j][t] = maxx
            else:
                delta[j][t] = maxx * B[j][int(seq[t])]

    viterbi_list1.append([0, '-', 0, delta[0][0]])

    while True:

        viterbi_list2 = []

        for item in range(len(viterbi_list1)):

            current_state = viterbi_list1[item][0]
            current_partial_seq = viterbi_list1[item][1]
            current_t = viterbi_list1[item][2]

            if current_state == 3*T+2 and current_t == T2-1:
                viterbi_list2.append(viterbi_list1[item])
                continue

            for j in range(3*T+3):
                if A[current_state][j] != 0:

                    if j%3 == 0 and current_t < T2:        # D states
                        viterbi_list2.append([j, current_partial_seq+'-', current_t, delta[j][current_t]])

                    elif j == 3*T+2 and current_t < T2-1:    # E state
                        viterbi_list2.append([j, current_partial_seq + '-', current_t+1, delta[j][current_t+1]])

                    elif j%3 == 2 and j != 3*T+2 and current_t < T2-1:      # M states
                        viterbi_list2.append(
                            [j, current_partial_seq + convert_NumToChar(seq_type, int(seq[current_t+1])),
                             current_t+1, delta[j][current_t+1]])

                    elif j%3 == 1 and current_t < T2-1:     # I states
                        viterbi_list2.append(
                            [j, current_partial_seq, current_t+1, delta[j][current_t+1]])

        viterbi_list1 = viterbi_list2
        viterbi_list1.sort(reverse=True, key=get_delta)
        viterbi_list1 = viterbi_list1[0:500]

        # print(viterbi_list1)

        terminate = True
        for item in range(len(viterbi_list1)):
            if int(viterbi_list1[item][2]) < T2-1:
                terminate = False
                break
        if terminate:
            break

    output = ""
    for item in range(len(viterbi_list1)):
        if int(viterbi_list1[item][0] == 3*T+2):
            output = viterbi_list1[item][1]
            break

    return output

# --------------------------------------------------------------------
#MAIN PROGRAM

addr = ""

if __name__ == "__main__":
    if len(sys.argv) < 2:
        addr = input('Enter dir of test cases: ')
    else:
        addr = sys.argv[1]

addr_in = addr + "/in"
addr_out = addr + "/out"

file_count = len([name for name in os.listdir(addr_in)
                  if os.path.isfile(os.path.join(addr_in, name))])

for file in range(file_count):

    f = open(addr_in+"/input"+str(file+1)+".txt", "r")
    list=[f.readline()[0:5]]
    for i in range(51):
        list.append(f.readline()[0:100])
    n, theta = list[0].split(' ')
    n = int(n)  # the number of seqs in MSA
    theta = int(theta)  # threshold for number of gaps in each column
    seq = list[1]
    old_seq_len = len(seq)-1
    old_seqs = np.zeros(shape=(n, old_seq_len), dtype='<U1')
    for j in range(old_seq_len):
        old_seqs[0][j] = seq[j]
    k=2
    for i in range(1, n):
        seq = list[k]
        for j in range(old_seq_len):
            old_seqs[i][j] = seq[j]
        k+=1
    seq = list[k-1]
    eval_seq = np.zeros(len(seq)-1, dtype='<U1')
    for j in range(len(eval_seq)):
        eval_seq[j] = seq[j]
    filtered_seqs = filtering(old_seqs, old_seq_len, theta)
    T = len(filtered_seqs[0])  # length of each seq after filtering
    seq_type = DNA_or_Prot(filtered_seqs, n, T)

    if seq_type:
        num_of_observations = 5
        # print('DNA!')
    else:
        num_of_observations = 21
        # print('Protein!')
    filtered_seqs2 = np.zeros(shape=(n, T+2), dtype=int)
    for i in range(n):
        for j in range(1, T):
            filtered_seqs2[i][j] = convert_charToNum(seq_type, filtered_seqs[i][j - 1])
        filtered_seqs2[i][0] = convert_charToNum(seq_type, '-')
        filtered_seqs2[i][T + 1] = convert_charToNum(seq_type, '-')
    print(len(filtered_seqs2))
    T = T + 2
    eval_seq2 = np.zeros(len(eval_seq)+2, dtype=int)
    for j in range(1, len(eval_seq)):
        eval_seq2[j] = convert_charToNum(seq_type, eval_seq[j - 1])
    eval_seq2[0] = convert_charToNum(seq_type, '-')
    eval_seq2[len(eval_seq) + 1] = convert_charToNum(seq_type, '-')

    T2 = len(eval_seq2)

    A = np.zeros(shape=(3*T+3, 3*T+3), dtype=float)
    B = np.zeros(shape=(3*T+3, num_of_observations), dtype=float)
    pi = np.zeros(3*T+3, dtype=float)
    A, B, pi = initialization(filtered_seqs2, A, B, pi, T, n, num_of_observations)
    print("--------------------Initialization Matrices-----------------")
    print("A: ",A,"\n")
    print("B: ",B,"\n")
    print("pi: ",pi,"\n")
    A, B, pi = Baum_Welch(filtered_seqs2, A, B, pi, T, n, num_of_observations)
    print("-------------Matrices After Baum Welch Function--------------")
    print("A: ",A,"\n")
    print("B: ",B,"\n")
    print("pi: ",pi,"\n")
    output = Viterbi2(eval_seq2, A, B, pi, T, T2, seq_type)
    output_seq = np.zeros(len(output) - 2, dtype='<U1')
    for i in range(0, len(output_seq)):
        output_seq[i] = output[i + 1]
    print("Output Sequence: ",output_seq,"\n")
    f2 = open(addr_out + "/output" + str(file+1) + ".txt", "r")
    right_output = f2.readline()
    right_output = right_output[0:len(right_output)-1]
    print("Right Output: ",right_output)
    str_of_input_seq = "".join(["".join(item) for item in eval_seq])
    str_of_output_seq = "".join(["".join(item) for item in output_seq])

    if str_of_output_seq == right_output:
        print(f"{str_of_input_seq} {right_output} {str_of_output_seq} True")
    else:
        print(f"{str_of_input_seq} {right_output} {str_of_output_seq} False")
