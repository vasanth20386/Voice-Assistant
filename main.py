class AminoAcid:  # DEFINING A CLASS FOR AMINO ACIDS, p_alpha & p_beta DENOTE THE PROPENSITY VALUES
    def __init__(self, name, abbr, p_alpha, p_beta):
        self.name = name
        self.abbr = abbr
        self.p_alpha = p_alpha
        self.p_beta = p_beta


a1 = AminoAcid('Alanine', 'A', 1.45, 0.97)  # NOW DEFINING ALL THE AMINO ACIDS
a2 = AminoAcid('Cysteine', 'C', 0.77, 1.30)
a3 = AminoAcid('Aspartic acid', 'D', 0.98, 0.80)
a4 = AminoAcid('Glutamic acid', 'E', 1.53, 0.26)
a5 = AminoAcid('Phenyl alanine', 'F', 1.12, 1.28)
a6 = AminoAcid('Glycine', 'G', 0.53, 0.81)
a7 = AminoAcid('Histidine', 'H', 1.24, 0.71)
a8 = AminoAcid('Isoleucine', 'I', 1, 1.6)
a9 = AminoAcid('Lysine', 'K', 1.07, 0.74)
a10 = AminoAcid('Leucine', 'L', 1.34, 1.22)
a11 = AminoAcid('Methionine', 'M', 1.2, 1.67)
a12 = AminoAcid('Asparagine', 'N', 0.73, 0.65)
a13 = AminoAcid('Proline', 'P', 0.59, 0.62)
a14 = AminoAcid('Glutamine', 'Q', 1.17, 1.23)
a15 = AminoAcid('Arginine', 'R', 0.79, 0.9)
a16 = AminoAcid('Serine', 'S', 0.79, 0.72)
a17 = AminoAcid('Threonine', 'T', 0.82, 1.2)
a18 = AminoAcid('Valine', 'V', 1.14, 1.65)
a19 = AminoAcid('Tryptophan', 'W', 1.14, 1.19)
a20 = AminoAcid('Tyrosine', 'Y', 0.61, 1.29)


def get_alpha_score(seq):  # FUNCTION THAT RETURNS ALPHA SCORE OF A SEQUENCE OF AMINO ACIDS
    alpha_score = 0
    i = 0
    while(i < len(seq)):
        if seq[i] == 'A':
            alpha_score += a1.p_alpha
        elif seq[i] == 'C':
            alpha_score += a2.p_alpha
        elif seq[i] == 'D':
            alpha_score += a3.p_alpha
        elif seq[i] == 'E':
            alpha_score += a4.p_alpha
        elif seq[i] == 'F':
            alpha_score += a5.p_alpha
        elif seq[i] == 'G':
            alpha_score += a6.p_alpha
        elif seq[i] == 'H':
            alpha_score += a7.p_alpha
        elif seq[i] == 'I':
            alpha_score += a8.p_alpha
        elif seq[i] == 'K':
            alpha_score += a9.p_alpha
        elif seq[i] == 'L':
            alpha_score += a10.p_alpha
        elif seq[i] == 'M':
            alpha_score += a11.p_alpha
        elif seq[i] == 'N':
            alpha_score += a12.p_alpha
        elif seq[i] == 'P':
            alpha_score += a13.p_alpha
        elif seq[i] == 'Q':
            alpha_score += a14.p_alpha
        elif seq[i] == 'R':
            alpha_score += a15.p_alpha
        elif seq[i] == 'S':
            alpha_score += a16.p_alpha
        elif seq[i] == 'T':
            alpha_score += a17.p_alpha
        elif seq[i] == 'V':
            alpha_score += a18.p_alpha
        elif seq[i] == 'W':
            alpha_score += a19.p_alpha
        else:
            alpha_score += a20.p_alpha
        i += 1
    return alpha_score


def get_beta_score(seq):  # FUNCTION THAT RETURNS BETA SCORE OF A SEQUENCE OF AMINO ACIDS
    beta_score = 0
    i = 0

    while (i < len(seq)):
        if seq[i] == 'A':
            beta_score += a1.p_beta
        elif seq[i] == 'C':
            beta_score += a2.p_beta
        elif seq[i] == 'D':
            beta_score += a3.p_beta
        elif seq[i] == 'E':
            beta_score += a4.p_beta
        elif seq[i] == 'F':
            beta_score += a5.p_beta
        elif seq[i] == 'G':
            beta_score += a6.p_beta
        elif seq[i] == 'H':
            beta_score += a7.p_beta
        elif seq[i] == 'I':
            beta_score += a8.p_beta
        elif seq[i] == 'K':
            beta_score += a9.p_beta
        elif seq[i] == 'L':
            beta_score += a10.p_beta
        elif seq[i] == 'M':
            beta_score += a11.p_beta
        elif seq[i] == 'N':
            beta_score += a12.p_beta
        elif seq[i] == 'P':
            beta_score += a13.p_beta
        elif seq[i] == 'Q':
            beta_score += a14.p_beta
        elif seq[i] == 'R':
            beta_score += a15.p_beta
        elif seq[i] == 'S':
            beta_score += a16.p_beta
        elif seq[i] == 'T':
            beta_score += a17.p_beta
        elif seq[i] == 'V':
            beta_score += a18.p_beta
        elif seq[i] == 'W':
            beta_score += a19.p_beta
        else:
            beta_score += a20.p_beta

        i+= 1

    return beta_score


def get_alpha_helix(seq):  # FUNCTION THAT RETURNS THE SITES WHERE ALPHA HELIX COULD BE FORMED
    s = []
    for i in range(len(seq)):
        s.append('_')
    if len(seq) < 6:
        print("Sequence length must be atleast 6 for this algorithm to work")
        return
    for i in range(len(seq)-5):
        n = 0
        for j in range(i, i+6):
            if get_alpha_score(seq[j]) >= 1:
                n += 1
        if n >= 4:
            for j in range(i, i + 6):
                if s[j] != 'H':
                    s[j] = 'H'
            p1 = i+6
            p2 = i-1
            while p1 < len(seq):
                if get_alpha_score(seq[p1-3:p1+1]) >= 4:
                    s[p1] = 'H'
                else:
                    break
                p1 += 1
            while p2 >= 0:
                if get_alpha_score(seq[p2:p2+4]) >= 4:
                    s[p2] = 'H'
                else:
                    break
                p2 -= 1
    alpha_helix = ""
    for i in range(len(s)):
        alpha_helix += s[i]
    return alpha_helix


def get_beta_sheet(seq):  # FUNCTION THAT RETURNS THE SITES WHERE BETA SHEETS COULD BE FORMED
    s = []
    for i in range(len(seq)):
        s.append('_')
    if len(seq) < 5:
        print("Sequence length must be atleast 5 for this algorithm to work")
        return
    for i in range(len(seq)-4):
        n = 0
        for j in range(i, i+5):
            if get_beta_score(seq[j]) >= 1:
                n += 1
        if n >= 3:
            for j in range(i, i+5):
                if s[j] != 'S':
                    s[j] = 'S'
            p1 = i+5
            p2 = i-1
            while p1 < len(seq):
                if get_beta_score(seq[p1-3:p1+1]) >= 4:
                    s[p1] = 'S'
                else:
                    break
                p1 += 1
            while p2 >= 0:
                if get_beta_score(seq[p2:p2+4]) >= 4:
                    s[p2] = 'S'
                else:
                    break
                p2 -= 1
    beta_sheet = ""
    for i in range(len(s)):
        beta_sheet += s[i]
    return beta_sheet


def conf_reso(seq1, seq2, seq):  # FUNCTION TO RESOLVE THE CONFLICT AT THE SITES WHERE BOTH ALPHA & BETA ARE FORMED
    res = []
    for i in range(len(seq)):
        res.append("")
    i = 0
    while i < len(seq):
        if (seq1[i] == 'H' and seq2[i] == '_') or (seq1[i] == '_' and seq2[i] == 'H'):
            res[i] = 'H'
            i += 1
        elif (seq1[i] == 'S' and seq2[i] == '_') or (seq1[i] == '_' and seq2[i] == 'S'):
            res[i] = 'S'
            i += 1
        elif seq1[i] == '_' and seq2[i] == '_':
            res[i] = 'T'
            i += 1
        else:
            n = 0
            while (seq1[i] == 'S' and seq2[i] == 'H') or (seq1[i] == 'H' and seq2[i] == 'S'):
                n += 1
                if i < len(seq)-1:
                    i += 1
                else:
                    break
            if i == len(seq)-1:
                i += 1
            p1 = get_alpha_score(seq[i-n:i])
            p2 = get_beta_score(seq[i-n:i])
            if p1 > p2:
                for k in range(i-n,i):
                    res[k] = 'H'
            else:
                for k in range(i-n,i):
                    res[k] = 'S'
    ans = ""
    for i in range(len(res)):
        ans += res[i]
    return ans


seq = "WHGCITVYWMTV" # INPUT SEQUENCE
seq1 = get_alpha_helix(seq)  # SITES WHERE ALPHA HELIXES ARE FORMED

seq2 = get_beta_sheet(seq)  # SITES WHERE BETA SHEETS ARE FORMED
final_seq = conf_reso(seq1, seq2, seq)  # FINAL SEQUENCE OBTAINED AFTER CONFLICT RESOLUTION

print("Input sequence is-")
print(seq)
print()
print("Sites where alpha-helixes are formed-")
print(seq1)
print()
print("Sites where beta sheets are formed-")
print(seq2)
print()
print("Final sequence obtained after conflict resolution-")
print(final_seq)