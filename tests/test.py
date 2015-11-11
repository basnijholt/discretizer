
import numpy as np
from collections import namedtuple

nt_test = namedtuple('nt_test', 'A P a')

par = nt_test(A = lambda x: 1.0 + (x > 1.5),
              P = lambda x: 1.0 + (x > 1.5),
              a = 1.0)


def hop_nn(xs, xt, par):
    A = par.A
    res = -(A(xt) + A(xs)) / (2*par.a**2)
    return res


def onsite(x, par):
    A = par.A
    res = (A(x-1)/2 + A(x) + A(x+1)/2) / (par.a**2)


sigma_0 = np.array([[1, 0], [0, 1]])
sigma_x = np.array([[0, 1], [1, 0]])
sigma_y = np.array([[0, -1j], [1j, 0]])
sigma_z = np.array([[1, 0], [0, -1]])


def hop_nn(xs, xt, par):
    A = par.A
    res = -((A(xt) + A(xs)) / (2*par.a**2) * sigma_0
            + np.sign(xs-xt) * 1j * sigma_x * (A(xs) + A(xt)) / (4 * par.a)
            - np.sign(xs-xt) * sigma_y * (A(xs) - A(xt)) / (4 * par.a))
    return res


def onsite(x, par):
    A = par.A
    res = (A(x-1)/2 + A(x) + A(x+1)/2) / (par.a**2) * sigma_0
    return res


par = nt_test(A = lambda x, y: 1.0 + (x > 1.5) + (y > 1.5),
              P = lambda x, y: 1.0 + (x > 1.5) + (y > 1.5),
              a = 1.0)


def hop_x_or_y(rs, rt, par):
    xs = rs[0]
    xt = rt[0]
    ys = rs[1]
    yt = rt[1]
    A = par.A
    res = -(A(xs, ys) + A(xt, yt)) / (2*par.a**2) * sigma_0
    return res


def hop_x_and_y(rs, rt, par):
    xs = rs[0]
    xt = rt[0]
    ys = rs[1]
    yt = rt[1]
    P = par.P
    res =  np.sign(xs-xt) * np.sign(ys-yt) * (P(xs, ys) + P(xt, yt)) * sigma_x / (4 * par.a**2)
    return res


def hop_x_or_y(rs, rt, par):
    xs = rs[0]
    xt = rt[0]
    ys = rs[1]
    yt = rt[1]
    A = par.A
    res = -(A(xs, ys) + A(xt, yt)) / (2*par.a**2) * sigma_0
    return res


def hop_x_and_y(rs, rt, par):
    xs = rs[0]
    xt = rt[0]
    ys = rs[1]
    yt = rt[1]
    P = par.P
