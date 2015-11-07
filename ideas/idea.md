# This is mostly to test how md works

## Simple case

$
\begin{aligned}
\hat{H} & = \frac{\hbar^2 k^2}{2m} \\
\hat{H} & = -\frac{\hbar^2}{2m} \partial_x^2
\end{aligned}
$

## Which in space in-dependent version gives

$
\hat{H} \varphi(x) = -\frac{\hbar^2}{2ma^2} \left[ \varphi(x+a) + \varphi(x-a) - 2\varphi(x) \right]
$

## let's try make it space dependent:

$
\begin{aligned}
\hat{H} & = k_x A(x) k_x \\
\hat{H} & = - \partial_x A(x) \partial_x \\
\end{aligned}
$

with
$ A(x)  = \frac{\hbar^2}{2m} $

### in one way approach it gives

$
\hat{H} \varphi(x) = -\frac{1}{a^2} \Big[ A(x+\frac{a}{2}) \varphi(x+a) + A(x-\frac{a}{2}) \varphi(x-a) - \big(A(x+\frac{a}{2}) + A(x-\frac{a}{2}) \big) \varphi(x) \Big]
$

hoppings:

$
\begin{aligned}
t_{x,\;x+a} & = -\frac{1}{a^2} A\big(x+\frac{a}{2}\big) \\
t_{x,\; x-a} & = -\frac{1}{a^2} A\big(x-\frac{a}{2}\big) \\
t_{x,\;x} & = -\frac{1}{a^2} \Big(A\big(x+\frac{a}{2}\big) + A\big(x-\frac{a}{2}\big)\Big) \\
\end{aligned}
$

### once in other (just hoppings)

hoppings:

$
\begin{aligned}
t_{x,\;x+a} & = -\frac{1}{a^2} \big[ + A(x+a) - A(x-a) + A(x)\big]  \\
t_{x,\;x-a} & = -\frac{1}{a^2} \big[ - A(x+a) + A(x-a) + A(x)\big]  \\
t_{x,\;x} & = -\frac{2}{a^2} A(x) \\
\end{aligned}
$
