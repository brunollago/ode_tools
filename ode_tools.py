import numpy as np
from sympy import *
from IPython.display import display, Math, Latex

class ode1:
    @staticmethod
    def solve_integrate_factor(a_1, a_0, g):
        x = symbols('x')
        y = Function('y')
        
        #step 1
        P_x = a_0/a_1
        f_x = g/a_1
        eq1 = Eq(y(x).diff(x) + P_x * y(x), f_x)      
        print("Escrevendo na forma padrão, temos:")
        display(eq1)

        #step 2
        FI_s = exp(Integral(P_x))
        FI = exp(integrate(P_x))
        FI_eq = Eq(FI_s, FI)
        print("O fator integrante é dado por:")
        display(FI_eq)

        #step 3
        step3 = Eq(FI * y(x).diff(x) + FI * x * y(x), FI * f_x)
        print("Multiplicando a equação pelo fator integrante, chegamos a:")
        display(Math(latex(step3)))

        #step 4
        step4 = Eq((y(x) * FI).diff(x, evaluate=False), FI * f_x)
        print("Identificando o lado esquerdo da equação como a derivada de um produto, podemos escrever:")
        display(Math(latex(step4)))

        #step 5
        step5 = Eq(y(x) * FI, Integral(FI * f_x))
        print("Agora podemos integrar a equação:")
        display(step5)
     
        #step 6
        c = symbols('c')
        step6 = Eq(y(x), (integrate(FI * f_x) + c) / FI)
        print("Por fim chegamos ao resultado:")
        display(step6)
    
    @staticmethod
    def solve_separable_variables(fr_x, gr_y, fl_x = 1, gl_y = 1):
        x = symbols('x')
        y = symbols('y')
        yy = Function('y')

        #step 0
        print("Podemos rescrever a equação da forma:")
        eq0 = Eq(yy(x).diff(x), fr_x * gr_y / fl_x / gl_y)
        display(eq0)

        #step 1
        dx, dy = symbols('dx dy')
        eq1 = Eq(fl_x * gl_y * dy, fr_x * gr_y * dx)
        print("Ou ainda:")
        display(eq1)

        #step 2
        step2 = Eq(gl_y / gr_y * dy, fr_x / fl_x * dx)
        print("Separando as variáveis, temos que:")
        display(step2)

        #step 3
        step3 = Eq(Integral(gl_y/gr_y), Integral(fr_x/fl_x))
        print("O próximo passo é integrar a equação:")
        display(step3)

        #step 4
        c = symbols('c')
        step4 = Eq(integrate(gl_y/gr_y), integrate(fr_x/fl_x) + c)
        print("O que resulta em:")
        display(step4)

        #step 5
        print("Tentamos encontrar uma expressão explícita de y(x):")
        step5 = Eq(y, solve(step4, y)[0])
        display(step5)

        #step 6
        print("E simplificar o resultado:")
        step6 = Eq(y, factor(solve(step4, y)[0]))
        display(step6)
    
    @staticmethod
    def solve_exact_equations(M_xy, N_xy):
        x = symbols('x')
        y = symbols('y')
        dx, dy = symbols('dx dy')
        
        #step 0
        eq0 = Eq(M_xy * dx + N_xy * dy, 0)
        print("Escrevendo na forma padrão, temos:")
        display(eq0)

        #step1
        ff = symbols('f')
        eq1 = Eq(ff.diff(x, evaluate=False), M_xy)
        print("Vamos começar utilizando M:")
        text = r"\dfrac{\partial f}{\partial x}=M(x,y)\rightarrow"
        display(Math(text + latex(eq1)))

        #step 2
        cc = symbols('c(y)')
        step2 = Eq(ff, integrate(M_xy, x) + cc)
        print("Integrando a equação, chegamos ao seguinte resultado:")
        display(step2)

        #step 3
        M_y = integrate(M_xy, x)
        step3 = Eq(cc.diff('y', evaluate=False), N_xy - M_y.diff(y))
        print("Agora derivamos em relação a y e utilizamos N:")
        text = r"\dfrac{d c}{d y}=N(x,y) -" + latex(M_y.diff(y))
        text = text + r"\rightarrow"
        display(Math(text + latex(step3)))

        #step 4
        c2 = symbols('c2')
        step4 = Eq(cc, integrate(N_xy - M_y.diff(y), y))
        print("Após a integração, encontramos o fator c:")
        display(step4)

        #step 5
        c = symbols('c')
        step5 = Eq(ff, integrate(M_xy, x) + integrate(N_xy - M_y.diff(y), y))
        print("E o resultado final é dado por:")
        display(step5)

class ode2:
    @staticmethod
    def solve_lin_h_coef_const(a, b, c):
        m = symbols('m')
        x = symbols('x')
        c1, c2 = symbols('c_1 c_2')

        #step 0
        print("A equação diferencial é a seguinte:")
        y2, y1, y = symbols('y\'\' y\' y')
        step0 = Eq(a*y2 + b*y1 + c*y, 0)
        display(step0)

        #step 1
        print("Neste caso, a equação auxiliar é:")
        step1 = Eq(a*m**2 + b*m + c, 0)
        display(step1)

        #step 2
        print("Com solução para m:")
        res = solve(step1, m)
        display(solve(step1, m))

        #step 3
        if len(res) == 1:
            step3 = Eq(y, c1*exp(res[0]*x) + c2*x*exp(res[0]*x))
            print("Neste caso, as raízes são iguais e temos:")
            display(step3)            
        elif(res[0].is_real):
            step3 = Eq(y, c1*exp(res[0]*x) + c2*exp(res[1]*x))
            print("Neste caso, as raízes são reais e distintas, logo:")
            display(step3)
        else:
            rr = re(res[0])
            ii = im(res[0])
            step3 = Eq(y, exp(rr * x) * (c1 * cos(ii * x) + c2 * cos(ii * x)))
            print("Neste caso, as raízes são complexas conjugadas, resultando em:")
            display(step3)