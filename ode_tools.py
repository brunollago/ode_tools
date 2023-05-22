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
        text = text + r"\rightarrow "
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
            step3 = Eq(y, exp(rr * x) * (c1 * sin(ii * x) + c2 * cos(ii * x)))
            print("Neste caso, as raízes são complexas conjugadas, resultando em:")
            display(step3)
    
    @staticmethod
    def solve_reduction_of_order(a_2, a_1, a_0, y_1):
        x = symbols('x')
        y = Function('y')

        #step 0
        P_x = a_1/a_2
        Q_x = a_0/a_2
        step0 = Eq(y(x).diff(x, x) + P_x * y(x).diff(x) + Q_x * y(x), 0)      
        print("Escrevendo na forma padrão, temos:")
        display(step0)

        #step 1
        text = r"Data a solução: " + latex(y_1)
        print("Dada a solução:")
        display(Math(r'y_1=' + latex(y_1)))

        #step 2
        print("Uma segunda solução L.I. é dada por:")
        sy_1, sy_2 = symbols('y_1 y_2')
        sdx, sppx = symbols('dx P(x)')
        left = sy_2
        middle = sy_1*Integral(exp(-Integral(sppx, x))/sy_1**2, x)
        right = y_1 * Integral(exp(-Integral(P_x,x))/y_1**2, x)
        text = latex(left) + "=" + latex(middle) +"="+ latex(right)
        display(Math(text))

        #step 3
        print("Resolvendo a integral chegamos a:")
        y_2 = y_1 * integrate(exp(-integrate(P_x,x))/y_1**2, x)
        step3 = Eq(sy_2, y_2)
        display(step3)

        #step 4
        print("E a solução geral é dada por:")
        step3 = Eq(symbols('y'), symbols('c_1')*y_1 + symbols('c_2')*y_2)
        display(step3)
    
    @staticmethod
    def solve_lin_h_cauchy_euler(a, b, c):
        m = symbols('m')
        x = symbols('x')
        c1, c2 = symbols('c_1 c_2')

        #step 0
        print("A equação diferencial é a seguinte:")
        y2, y1, y = symbols('y\'\' y\' y')
        step0 = Eq(a*x**2*y2 + b*x*y1 + c*y, 0)
        display(step0)

        #step 1
        print("Neste caso, a equação auxiliar é:")
        step1 = Eq(a*m*(m-1) + b*m + c, 0)
        display(step1)

        #step 2
        print("Com solução para m:")
        res = solve(step1, m)
        display(solve(step1, m))

        #step 3
        if len(res) == 1:
            step3 = Eq(y, c1*x**res[0] + c2*ln(x)*x**res[0])
            print("Neste caso, as raízes são iguais e temos:")
            display(step3)            
        elif(res[0].is_real):
            step3 = Eq(y, c1*x**res[0] + c2*x**res[1])
            print("Neste caso, as raízes são reais e distintas, logo:")
            display(step3)
        else:
            rr = re(res[0])
            ii = im(res[0])
            step3 = Eq(y, x**rr * (c1 * cos(ii * ln(x)) + c2 * cos(ii * ln(x))))
            print("Neste caso, as raízes são complexas conjugadas, resultando em:")
            display(step3)
    
    @staticmethod
    def solve_lin_nh_variation_of_parameters(a_2, a_1, a_0, g_x, y_1, y_2):
        x = symbols('x')

        #step 0
        y2, y1, y = symbols('y\'\' y\' y')
        P_x = a_1 / a_2
        Q_x = a_0 / a_2
        f_x = g_x / a_2
        step0 = Eq(y2 + P_x * y1 + Q_x * y, f_x)
        print("A equação é da forma:")
        display(step0)

        #step 0.5
        print("Dadas as soluções da homogênea, a solução particular é da forma:")
        text = r"y_p=u_1y_1+u_2y_2"
        display(Math(text))
        print("Onde:")
        text = r"u'_1=\dfrac{\rm{det}(W_1)}{\rm{det}(W)}\quad ;\quad u'_2=\dfrac{\rm{det}(W_2)}{\rm{det}(W)}"
        display(Math(text))
        print("com:")
        W = Matrix([[y_1, y_2],[y_1.diff(x), y_2.diff(x)]])
        W_1 = Matrix([[0, y_2],[f_x, y_2.diff(x)]])
        W_2 = Matrix([[y_1, 0],[y_1.diff(x), f_x]])
        display(Math(r"W=" + latex(W)))
        display(Math(r"W_1=" + latex(W_1)))
        display(Math(r"W_2=" + latex(W_2)))

        #step 1 and 2
        step1_1 = Eq(symbols('u\'_1'), simplify(W_1.det() / W.det()))
        step1_2 = Eq(symbols('u\'_2'), simplify(W_2.det() / W.det()))
        step2_1 = Eq(symbols('u_1'), Integral(simplify(W_1.det() / W.det()), x))
        step2_2 = Eq(symbols('u_2'), Integral(simplify(W_2.det() / W.det()), x))
        print("As equações para u_1 e u_2 são:")
        display(Math(latex(step1_1) + r"\rightarrow " + latex(step2_1)))
        display(Math(latex(step1_2) + r"\rightarrow " + latex(step2_2)))

        #step 3
        u_1 = integrate(simplify(W_1.det() / W.det()), x)
        u_2 = integrate(simplify(W_2.det() / W.det()), x)
        step3_1 = Eq(symbols('u_1'), u_1)
        step3_2 = Eq(symbols('u_2'), u_2)
        print("Resolvendo as integrais, temos que:")
        display(step3_1)
        display(step3_2)

        #step 4
        c_1, c_2 = symbols('c_1 c_2')
        step4 = Eq(symbols('y'), simplify(c_1 * y_1 + c_2 * y_2 + u_1 * y_1 + u_2 * y_2)) # May cause problems related to constants
        print("E a solução geral é a seguinte:")
        display(step4)