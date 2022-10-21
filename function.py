import numpy as np

def SVD(S,f):
    U, sigma, Vt = np.linalg.svd(S)
    U = U[ : , : f ]
    sigma = np.diag(sigma[0:f])
    X = np.dot(U, sigma ** 0.5)
    return X

def cos_sim(vector1, vector2):
    dot_product = 0.0
    normA = 0.0
    normB = 0.0
    for a, b in zip(vector1, vector2):
        dot_product += a * b
        # print(dot_product)
        normA += a ** 2
        normB += b ** 2
    if normB == 0.0 or normA == 0.0:
        return 0
    else:
        return dot_product / ((normA ** 0.5) * (normB ** 0.5))

def lapsim(dd):
    mm = np.zeros(dd.shape)
    for i in range(mm.shape[0]):
        for j in range(mm.shape[1]):
            mm[i, j] = dd[i, j] / np.sqrt(np.sum(dd[i, :]) * np.sum(dd[j, :]))
    return mm

def tpr_fpr_precision_recall(true, pred):
    tp, fp, tn, fn = 0, 0, 0, 0
    index = list(reversed(np.argsort(pred)))
    tpr = []
    fpr = []
    precision = []
    recall = []
    for i in range(pred.shape[0]):
        if true[int(index[i])] == 1:
            tp += 1
        else:
            fp += 1
        if np.sum(true) == 0:
            tpr.append(0)
            fpr.append(0)
            precision.append(0)
            recall.append(0)
        else:
            tpr.append(tp / np.sum(true))
            fpr.append(fp / (true.shape[0] - np.sum(true)))
            precision.append(tp / (tp + fp))
            recall.append(tp / np.sum(true))
    return tpr, fpr, precision, recall

def equal_len_list(a):
    row_len = []
    for i in a:
        row_len.append(len(i))
    min_len = min(row_len)
    equal_len_a = []
    for i in a:
        tem_list = []
        multi = len(i)/min_len
        for j in range(min_len-1):
            tem_list.append(i[int(j*multi)])
        tem_list.append(i[-1])
        equal_len_a.append(tem_list)
    return equal_len_a

class IKNN:
    def __init__(self, k) -> None:
        self._k = k

    def fit(self, _X):
        self.X = _X.copy()
        return self

    def neighbors(self, i) -> np.ndarray:
        drugs_sim = self.X[i]
        values = np.sort(drugs_sim)[::-1][1:self._k + 1]
        indexes = np.argsort(drugs_sim)[::-1][1:self._k + 1]
        return (indexes, values)


def pruning(K, drug_mat, target_mat, miu):
    
    m = len(drug_mat)
    n = len(target_mat)
    sr= np.zeros((m, m))
    sp= np.zeros((n, n))

    iknn = IKNN(K)
    weights = np.zeros(K)

    iknn.fit(drug_mat)
    for d in range(m):
        (indexes, values) = iknn.neighbors(d)
        z = np.sum(values)
        if z==0:
            continue
        for i in range(K):
            sr[d,indexes[i]] = (miu ** i) * values[i]

    iknn.fit(target_mat)
    for t in range(n):
        (indexes, values) = iknn.neighbors(t)
        z = np.sum(values)
        if z==0:
            continue
        for i in range(K):
            sp[indexes[i], t] = (miu ** i) * values[i] 
    return sr, sp



def BiRW(R, sr, sp, a, rs):
    sr=drug_rowsimm(sr)
    sp=dis_colsimm(sp)
    R1 = rs * sr.dot(R) + (1 - rs) * a
    R2 = rs * R.dot(sp) + (1 - rs) * a
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            R[i][j] = max(R[i][j], (R1[i][j] + R2[i][j]) / 2)
    return R


def drug_rowsimm(dd):
    mm = np.zeros(dd.shape)
    for i in range(dd.shape[0]):
        # for j in range(mm.shape[1]):
            if np.sum(dd[i, :])!=0:
               mm[i, :] = dd[i, :] / np.sum(dd[i, :])
    return mm



def dis_colsimm(dd):
    mm = np.zeros(dd.shape)
    for j in range(dd.shape[1]):
            if np.sum(dd[:, j]) != 0:
                mm[:, j] = dd[:, j] / np.sum(dd[:, j])
    return mm

