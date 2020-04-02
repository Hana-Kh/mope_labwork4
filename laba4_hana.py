import random
import numpy
from scipy.stats import t,f


#create new y array
def new_y(n, m, ymin, ymax):
    y = [[random.randint(ymin, ymax) for i in range(m)] for k in range(n)]
    return y


#add y column
def append_y(n, y, ymin, ymax):
    for i in y:
        i.append(random.randint(ymin, ymax))


def my(ymat):
    ymmat = [sum(i)/3 for i in ymat]
    return ymmat


def get_b_norm(xnorm, ym):
    n = len(ym)
    b = [0 for i in range(n)]
    b[0] = sum(ym) / n
    for i in range(n):
        b[1] += ym[i] * xnorm[i][1] / n
        b[2] += ym[i] * xnorm[i][2] / n
        b[3] += ym[i] * xnorm[i][3] / n
        if n == 8:
            b[4] += ym[i] * xnorm[i][1] * xnorm[i][2] / n
            b[5] += ym[i] * xnorm[i][1] * xnorm[i][3] / n
            b[6] += ym[i] * xnorm[i][2] * xnorm[i][3] / n
            b[7] += ym[i] * xnorm[i][1] * xnorm[i][2] * xnorm[i][3] / n
    return b


def get_b_nat(m, xnat, ym):
    n = len(ym)
    def ai1(x, k):
        a = [0 for i in range(n)]
        for i in range(n):
            a[0] += x[i][k]
            a[1] += xnat[i][0] * x[i][k]
            a[2] += xnat[i][1] * x[i][k]
            a[3] += xnat[i][2] * x[i][k]
            if n == 8:
                a[4] += xnat[i][0] * xnat[i][1] * x[i][k]
                a[5] += xnat[i][0] * xnat[i][2] * x[i][k]
                a[6] += xnat[i][1] * xnat[i][2] * x[i][k]
                a[7] += xnat[i][0] * xnat[i][1] * xnat[i][2] * x[i][k]
        return a

    def ai2(x, k, l):
        a = [0 for i in range(n)]
        for i in range(n):
            a[0] += x[i][k] * x[i][l]
            a[1] += xnat[i][0] * x[i][k] * x[i][l]
            a[2] += xnat[i][1] * x[i][k] * x[i][l]
            a[3] += xnat[i][2] * x[i][k] * x[i][l]
            if n == 8:
                a[4] += xnat[i][0] * xnat[i][1] * x[i][k] * x[i][l]
                a[5] += xnat[i][0] * xnat[i][2] * x[i][k] * x[i][l]
                a[6] += xnat[i][1] * xnat[i][2] * x[i][k] * x[i][l]
                a[7] += xnat[i][0] * xnat[i][1] * xnat[i][2] * x[i][k] * x[i][l]
        return a

    def ai3(x, k, l, m):
        a = [0 for i in range(n)]
        for i in range(n):
            a[0] += x[i][k] * x[i][l] * x[i][m]
            a[1] += xnat[i][0] * x[i][k] * x[i][l] * x[i][m]
            a[2] += xnat[i][1] * x[i][k] * x[i][l] * x[i][m]
            a[3] += xnat[i][2] * x[i][k] * x[i][l] * x[i][m]
            if n == 8:
                a[4] += xnat[i][0] * xnat[i][1] * x[i][k] * x[i][l] * x[i][m]
                a[5] += xnat[i][0] * xnat[i][2] * x[i][k] * x[i][l] * x[i][m]
                a[6] += xnat[i][1] * xnat[i][2] * x[i][k] * x[i][l] * x[i][m]
                a[7] += xnat[i][0] * xnat[i][1] * xnat[i][2] * x[i][k] * x[i][l] * x[i][m]
        return a

    a = []
    a1 = [0 for i in range(n)]
    a1[0] = n
    for i in range(n):
        a1[1] += xnat[i][0]
        a1[2] += xnat[i][1]
        a1[3] += xnat[i][2]
        if n == 8:
            a1[4] += xnat[i][0] * xnat[i][1]
            a1[5] += xnat[i][0] * xnat[i][2]
            a1[6] += xnat[i][1] * xnat[i][2]
            a1[7] += xnat[i][0] * xnat[i][1] * xnat[i][2]

    a.append(a1)
    a.append(ai1(xnat, 0))
    a.append(ai1(xnat, 1))
    a.append(ai1(xnat, 2))
    if n == 8:
        a.append(ai2(xnat, 0, 1))
        a.append(ai2(xnat, 0, 2))
        a.append(ai2(xnat, 1, 2))
        a.append(ai3(xnat, 0, 1, 2))

    c = [0 for i in range(n)]
    for i in range(n):
        c[0] += ym[i]
        c[1] += ym[i] * xnat[i][0]
        c[2] += ym[i] * xnat[i][1]
        c[3] += ym[i] * xnat[i][2]
        if n == 8:
            c[4] += ym[i] * xnat[i][0] * xnat[i][1]
            c[5] += ym[i] * xnat[i][0] * xnat[i][2]
            c[6] += ym[i] * xnat[i][1] * xnat[i][2]
            c[7] += ym[i] * xnat[i][0] * xnat[i][1] * xnat[i][2]

    ax = numpy.array(a)
     
    cx = numpy.array(c)

    b = numpy.linalg.solve(ax, cx)
    return b


def table_student(prob, n, m):
    x_vec = [i*0.0001 for i in range(int(5/0.0001))]
    par = 0.5 + prob/0.1*0.05
    f3 = (m - 1) * n
    for i in x_vec:
        if abs(t.cdf(i, f3) - par) < 0.000005:
            return i


def table_fisher(prob, n, m, d):
    x_vec = [i*0.001 for i in range(int(10/0.001))]
    f3 = (m - 1) * n
    for i in x_vec:
        if abs(f.cdf(i, n-d, f3)-prob) < 0.0001:
            return i


def comb(arr):
    return [1, *arr, arr[0]*arr[1], arr[0]*arr[2], arr[1]*arr[2], arr[0]*arr[1]*arr[2]]


def kohren(n, m, prob, disp):
    fisher = table_fisher(prob, n, m, 1)
    gt = fisher/(fisher+(m-1)-2)
    return max(disp) / sum(disp) < gt


def student(m, prob, disp, xnorm, ym):
    n = len(ym)
    sbt = (sum(disp) / m) ** (0.5) / n
    
    beta = [sum([comb(xnorm[j][1:])[i] * ym[j] / n for j in range(n)]) for i in range(n)]
    
    tt = table_student(prob, n, m)
    st = [(abs(i) / sbt) > tt for i in beta]
    return st


def fisher(m, prob, disp, ym, xnat, b, d):
    n = len(ym)
    if d == n:
        return False
    
    sad = sum([(sum([comb(xnat[i])[j] * b[j] for j in range(n)]) - ym[i]) ** 2 for i in range(n)])        
    sad = sad * m / (n - d)
    fp = sad / sum(disp) / n
    ft = table_fisher(prob, n, m, d)
    
    return fp < ft


def console_print():
    titles_x = ["№", "X0", "X1", "X2", "X3", "X1*X2", "X1*X3", "X2*X3", "X1*X2*X3"]

    # cycle for pretty printing title of table with normal parameters
    for j in range(N+1):
            s = ""
            if j == 0:
                s = "  {:1s}  " # for №
            if j == 1:
                s = " {:2s}  "  # for X0
            if j >= 2 and j < 5:
                s = "  {}  "    # for X + num
            if j >= 5 and j < 8:
                s = " {:5s}  "  # for X*X, with different combinations
            if j == 8:
                s = " {:8s}  "  # for X*X*X
            print(s.format(titles_x[j]), end="")    # taking all titles from list

    # this cycle is used for printing Yi in title of table
    for i in range(m):
        print(" Yi{:d}  ".format(i+1), end="")

    # printing Y middle, Y experimental and dispersion
    print("   Ys       Ye       S^2   ", end="")

    print()
    # fill table with data
    for i in range(N):
        print("  {:1d}   {:2d}   {:3d}   {:3d}   {:3d}  ".format(i+1, xnorm[i][0], *xnat[i]), end="")
        if N == 8:
            print(" {:5d}   {:5d}   {:5d}   {:8d}  "
                  .format(xnat[i][0]*xnat[i][1], xnat[i][0]*xnat[i][2], xnat[i][1]*xnat[i][2],
                          xnat[i][0]*xnat[i][1]*xnat[i][2]), end="")
        for j in y[i]:
            print(" {:3d}  ".format(j), end="")
            
        yss = b[0]*d_arr[0] + b[1]*xnat[i][0]*d_arr[1] + b[2]*xnat[i][1]*d_arr[2] + b[3]*xnat[i][2]*d_arr[3]
        if N == 8:
            yss += b[4]*xnat[i][0]*xnat[i][1]*d_arr[4] + b[5]*xnat[i][0]*xnat[i][2]*d_arr[5] + b[6]*xnat[i][1]*xnat[i][2]*d_arr[6] + b[7]*xnat[i][0]*xnat[i][1]*xnat[i][2]*d_arr[7]
        
        print(" {:6.2f}   {:6.2f}   {:6.2f}  ".format(ym[i], yss, disp[i]), end="")
        
        print()

    print("\nNatural linear regrecy:  Y = ", end="")
    if d_arr[0] != 0:
        print("{:.5f}".format(b[0]), end="")
    for i in range(1, N):
        if d_arr[i] != 0:
            print(" + {:.5f}*{}".format(b[i], titles_x[i+1]), end="")

    print("\n")

    # the same style of printing table with natural parameters
    for j in range(N+1):
        s = ""
        if j == 0:
            s += "  {}  "       # for №
        if j == 1:
            s += " {}  "        # for X0
        if j >= 2 and j < 5:    
            s += " {}  "        # for X + num
        if j >= 5 and j < 8:
            s += " {}  "        # for X*X, with different combinations
        if j == 8:
            s += " {}  "        # for X*X*X
        print(s.format(titles_x[j]), end="")    # taking all titles from list

    # this cycle is used for printing Yi in title of table
    for i in range(m):
        print(" Yi{:d}  ".format(i+1), end="")

    # printing Y middle, Y experimental and dispersion
    print("   Ys       Ye       S^2   ", end="")

    print()
    # fill table with data
    for i in range(N):
        print("  {:1d}   {:2d}   {:2d}   {:2d}   {:2d}  ".format(i+1, *xnorm[i]), end="")
        if N == 8:
            print(" {:5d}   {:5d}   {:5d}   {:8d}  "
                  .format(xnorm[i][1]*xnorm[i][2], xnorm[i][1]*xnorm[i][3], xnorm[i][2]*xnorm[i][3],
                          xnorm[i][1]*xnorm[i][2]*xnorm[i][3]), end="")
        for j in y[i]:
            print(" {:3d}  ".format(j), end="")
            
        yss = bnorm[0]*d_arr[0] + bnorm[1]*xnorm[i][1]*d_arr[1] + bnorm[2]*xnorm[i][2]*d_arr[2] + bnorm[3]*xnorm[i][3]*d_arr[3]
        if N == 8:
            yss += bnorm[4]*xnorm[i][1]*xnorm[i][2]*d_arr[4] + bnorm[5]*xnorm[i][1]*xnorm[i][3]*d_arr[5] + bnorm[6]*xnorm[i][2]*xnorm[i][3]*d_arr[6] + bnorm[7]*xnorm[i][1]*xnorm[i][2]*xnorm[i][3]*d_arr[7]
        
        print(" {:6.2f}   {:6.2f}   {:6.2f}  ".format(ym[i], yss, disp[i]), end="")
        print()

    print("\nNormal linear regrecy:  Y = ", end="")
    if d_arr[0] != 0:
        print("{:.5f}".format(bnorm[0]), end="")
    for i in range(1, N):
        if d_arr[i] != 0:
            print(" + {:.5f}*{}".format(bnorm[i], titles_x[i+1]), end="")


m = 3
N = 4
prob = 0.95

x1min = -5
x1max = 15
x2min = -25
x2max = 10
x3min = -5
x3max = 20

xmmin = (x1min + x2min + x3min) / 3
xmmax = (x1max + x2max + x3max) / 3
ymax = round(200 + xmmax)
ymin = round(200 + xmmin)

xnorm = [[1, -1, -1, -1],
         [1, -1, 1, 1],
         [1, 1, -1, 1],
         [1, 1, 1, -1],
         [1, -1, -1, 1],
         [1, -1, 1, -1],
         [1, 1, -1, -1],
         [1, 1, 1, 1]]

xnat = [[x1min, x2min, x3min],
        [x1min, x2max, x3max],
        [x1max, x2min, x3max],
        [x1max, x2max, x3min],
        [x1min, x2min, x3max],
        [x1min, x2max, x3min],
        [x1max, x2min, x3min],
        [x1max, x2max, x3max]]

while True:
    print("\n\nStart N = {}".format(N))
    y = new_y(N, m, ymin, ymax)
    while True:
        print("\nm = {}\n".format(m))
        ym = my(y)

        b = get_b_nat(m, xnat, ym)
        bnorm = get_b_norm(xnorm, ym)

        disp = []
        for i in range(len(ym)):
            s = 0
            for k in range(m):
                s += (ym[i] - y[i][k]) ** 2
            disp.append(s / m)
        
        koh = kohren(N, m, prob, disp)
        print("Dispersion uniform is {}, with probability = {:.2}".format(koh, prob))
        if koh:
            break
        else:
            m += 1
            append_y(N, y, ymin, ymax)

    d_arr = student(m, prob, disp, xnorm, ym)
    d = sum(d_arr)
    
    fis = fisher(m, prob, disp, ym, xnat, b, d)
    print("Equation adequativity is {}, with probability = {:.2f}\n".format(fis, prob))
    
    console_print()
    if fis:
        break
    else:
        if N == 4:
            N = 8
        else:
            N = 4
        m = 3

        



