import numpy as np


def cleanData(data, diag_len=4):
    len_i = len(data[:, 0])
    len_j = len(data[0, :])

    new_data = np.zeros_like(data)

    max_list = np.argmax(data, axis=1)
    for i in range(len(max_list)):
        new_data[i][max_list[i]] = 1

    new_data = np.zeros_like(data)

    max_list = np.argmax(data, axis=0)
    for i in range(len(max_list)):
        if data[max_list[i]][i] > 0:
            new_data[:, i] = 0
            new_data[max_list[i], i] = 1

    diag_len_by_two = int(diag_len / 2)
    score_density = np.copy(new_data)
    score_copy = np.copy(score_density)

    for i in range(diag_len + 2, len_i - diag_len - 1):
        for j in range(diag_len + 2, len_j - diag_len - 1):
            value = 0
            for w in range(-diag_len_by_two, diag_len_by_two + 1):
                if score_density[i + w, j + w] > 0:
                    value = value + 1
            if value >= diag_len:
                for k in range(1, diag_len + 2):
                    score_copy[i + k, j + k] = 1
                    score_copy[i - k, j - k] = 1
            elif score_copy[i, j] == 0:
                score_copy[i, j] = 0

    for i in range(diag_len + 2, len_i - diag_len - 1):
        for j in range(diag_len + 2, len_i - diag_len - 1):
            value = 0
            for w in range(-diag_len_by_two, (diag_len_by_two + 1)):
                if score_density[i - w, j + w] > 0:
                    value = value + 1

            if value >= diag_len:
                for k in range(1, diag_len + 2):
                    score_copy[i - k, j + k] = 1
                    score_copy[i + k, j - k] = 1
            elif score_copy[i, j] == 0:
                score_copy[i, j] = 0

    for i in range(0, (len(score_copy[:, 1]))):
        for j in range(0, (len(score_copy[1, :]))):
            # value = 0
            min_i = max(0, i - 1)
            max_i = min(len_i, i + 2)
            min_j = max(0, j - 1)
            max_j = min(len_j, j + 2)
            value = np.sum(score_copy[min_i:max_i, min_j:max_j])
            if value < 2:
                score_copy[i, j] = 0
    return score_copy


def calculateScore(data):
    len_i = len(data[:, 0])
    len_j = len(data[0, :])

    score = 0
    dist_th = 1.5
    dvec1 = abs(np.argmax(data[:, 1]) - np.argmax(data[:, 0]))
    dvec2 = abs(np.argmax(data[:, 2]) - np.argmax(data[:, 1]))
    dvec3 = abs(np.argmax(data[:, 3]) - np.argmax(data[:, 2]))

    for i in range(4, len_i):
        distance = np.mean([dvec1, dvec2, dvec3])

        dvec1 = dvec2
        dvec2 = dvec3
        dvec3 = abs(np.argmax(data[:, i]) -
                    np.argmax(data[:, i - 1]))
        if (distance > dist_th or distance == 0):
            score = score + len_i

    score = (score / (len_i ** 2))
    print(score)
    return score


def downsample(matrix, downscale):
    mat = np.copy(matrix)
    l = min(len(mat[1, :]), len(mat[:, 1]))
    size = int(l/downscale)
    m = np.zeros((size, size), dtype=np.int)
    for i in range(1, l):
        for j in range(1, l):
            mup = max(1, i-1)
            mdown = min(l, i+1)
            mleft = max(1, j-1)
            mright = min(l, j+1)

            if(np.sum(mat[mup:mdown, mleft:mright]) > 0):
                m[int(i/downscale), int(j/downscale)] = 1
    return (m)


def growing_regions(mat, reward=6, penalty=15, sidePenalty=3, MAXHSPS=500, TH=5, WSIZE=7):
    l = len(mat[1, ])
    HSPS = np.zeros((MAXHSPS, 5), dtype=np.int)

    idx = 1
    lH = round(WSIZE/2) - 1
    rH = round(WSIZE/2) + 1
    if(WSIZE % 2 == 0):
        print("WSIZE MUST BE ODD")
        return
    i = 0
    while(i < l-1):
        value = max(mat[i, ]) * reward
        if(value == 0):
            i = i + 1
        pos = np.argmax(mat[i, ])
        endfrag = pos
        j = i
        count_penalties = 1
        while(value > 0 and j < l-1):
            mat[max(0, j-1), max(0, endfrag-2)] = 0
            mat[max(0, j-1), max(0, endfrag-1)] = 0
            mat[max(0, j-1), endfrag] = 0
            mat[max(0, j-1), min(l-1, endfrag+1)] = 0
            mat[max(0, j-1), min(l-1, endfrag+2)] = 0
            mat[j, max(0, endfrag-2)] = 0
            mat[j, max(0, endfrag-1)] = 0
            mat[j, endfrag] = 0
            mat[j, min(l-1, endfrag+1)] = 0
            mat[j, min(l-1, endfrag+2)] = 0
            j = j + 1
            mleft = max(0, endfrag-lH)
            mright = min(l-1, endfrag+lH+1)
            window = mat[j, mleft:mright]
            v = max(window)
            selected = np.argmax(window)
            chose_diagonal = False
            if(len(window) == WSIZE and v == window[lH]):
                selected = lH
                chose_diagonal = True
            if(len(window) == WSIZE and v == window[rH]):
                selected = rH
                chose_diagonal = True
            if(v != 0):
                endfrag = (mleft + selected)  # To make the indexing
                if(len(window) == WSIZE):
                    endfrag = endfrag - 1
                endfrag = max(1, min(l, endfrag))
            if(v == 0):
                value = value - count_penalties * penalty
                count_penalties = count_penalties + 1
            else:
                if(not chose_diagonal):
                    value = value + count_penalties * (-sidePenalty)
                    count_penalties = count_penalties + 1
                else:
                    count_penalties = 1
                    value = value + reward
        if(j-i > TH):
            HSPS[idx, 0] = pos
            HSPS[idx, 1] = i
            HSPS[idx, 2] = endfrag
            HSPS[idx, 3] = j
            HSPS[idx, 4] = abs(i-j)
            idx = idx + 1
        if(idx == MAXHSPS):
            break
    return (HSPS)


def detect_events(HSPS, sampling):
    DIAG_SEPARATION = 10
    # same as HSPS but adding the event
    output = np.zeros((len(HSPS[:, 1]), 1+len(HSPS[1, :])))
    event_types = ['' for _ in range(len(HSPS))]

    j = 0
    for i in range((len(HSPS[:, 1]))):
        if(sum(HSPS[i, :]) > 0):
            j = j + 1
            is_inverted = False
            is_diagonal = True

            if(HSPS[i, 0] > HSPS[i, 2]):
                is_inverted = True
            if(abs(HSPS[i, 0] - HSPS[i, 1]) > DIAG_SEPARATION and abs(HSPS[i, 2] - HSPS[i, 3]) > DIAG_SEPARATION):
                is_diagonal = False

            output[i, 0] = HSPS[i, 0] * sampling
            output[i, 1] = HSPS[i, 1] * sampling
            output[i, 2] = HSPS[i, 2] * sampling
            output[i, 3] = HSPS[i, 3] * sampling
            output[i, 4] = HSPS[i, 4] * sampling

            if(is_diagonal):
                event_types[i] = "synteny block"
            if(is_diagonal and is_inverted):
                event_types[i] = "inversion"
            if(not is_diagonal and not is_inverted):
                event_types[i] = "transposition"
            if(not is_diagonal and is_inverted):
                event_types[i] = "inverted transposition"
    result = []
    for i in range(1, j+1):
        row = []
        for k in range(5):
            row.append(output[i, k])
        row.append(event_types[i])
        result.append(row)
    return result
