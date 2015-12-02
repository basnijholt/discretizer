def follow_path(expr, path):
    res = expr
    for i in np.arange(len(path)):
        res = res.args[path[i]]
    return res

def interchange(expr, sub, path):
    res = sub
    for i in np.arange(len(path)):
        temp = follow_path(expr, path[:-(i+1)])
        args = list(temp.args)
        args[path[len(path)-i-1]] = res
        res = temp.func(*tuple(args))
    return res

def interpolate_Function(expr):
    path = 'None'
    factor = 'None'
    change = False
    summand_0 = 'None'
    summand_1 = 'None'
    res = expr
    for i in np.arange(len(expr.args)):
        argument = sympy.expand(expr.args[i])
        if argument.func == sympy.Add:
            for j in np.arange(len(argument.args)):
                summand = argument.args[j]
                if summand.func == sympy.Mul:
                    for k in np.arange(len(summand.args)):
                        temp = 0
                        if summand.args[k] == sympy.Symbol('a'):
                            temp = sympy.Mul(sympy.Mul(*summand.args[:k]),
                                             sympy.Mul(*summand.args[k+1:]))
                            #print(temp)
                        if not temp == int(temp):
                            #print('found one')
                            factor = (temp)
                            path = np.array([i, j, k])
    if not factor == 'None':
        change = True

        sign = np.sign(factor)
        offsets = np.array([int(factor), sign * (int(sign * factor) + 1)])
        weights = 1/np.abs(offsets - factor)
        weights = weights/np.sum(weights)

        res = (  weights[0] * interchange(expr, offsets[0] * sympy.Symbol('a'), path[:-1])
               + weights[1] * interchange(expr, offsets[1] * sympy.Symbol('a'), path[:-1]))

    return sympy.expand(res), change

def interpolate(expr):
    change = False
    expr = sympy.expand(expr)
    res = expr

    if isinstance(expr, sympy.Function):# and not change:
        path = np.array([])
        temp, change = interpolate_Function(follow_path(expr, path))
        res = interchange(expr, temp, path)

    for i in np.arange(len(expr.args)):
        path = np.array([i])
        if isinstance(follow_path(expr, path), sympy.Function) and not change:
            temp, change = interpolate_Function(follow_path(expr, path))
            res = interchange(expr, temp, path)

        for j in np.arange(len(expr.args[i].args)):
            path = np.array([i, j])
            if isinstance(follow_path(expr, path), sympy.Function) and not change:
                temp, change = interpolate_Function(follow_path(expr, path))
                res = interchange(expr, temp, path)

    if change:
        res = interpolate(res)

    return sympy.expand(res)