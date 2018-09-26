# Wissenschaftliche Simulation Projekt
# Projektaufgabe 7 : Integralgleichung
# Date:  19.07.2017
# Autor: Abdellah Frindou    10006874
#        Abdessatar ben ifa  10007370
#        Ayoub Bahaj         10006875



# model_F : Function calculate f(x)
# and Stock it Vector
function model_F(x::Vector{Float64})
        f= zeros(n)
        for i=1:length(f)
            sum = 0
                for j=1:n
                    sum = sum+(cos((i-0.5)*(j-0.5)*(h^2))*(x[j]^3))
                end
          f[i] = x[i] + h * sum - 2;
        end
        return f
end

#model_DF : Function Calculate and fill Jacobi-Matrix
function model_DF(x::Vector{Float64})
        A= zeros(n,n)
          for i=1:n
              for j=1:n

                    if ( i==j)
                        A[i,j] = 1 + 1/20*cos((i-1/2)^2*(h^2))*(x[j]^2.0)
                    else
                        A[i,j] = 1/20*cos((i-1/2)*(j-1/2)*(h^2))*(x[j]^2.0)
                    end
              end
          end
        return A
end

# Gauss Elemination von normaler newton verfahren and Calculate SK
# the parameter bol is just a boolean to distinguish between the normaler Newton-Verfahren and vereinfachte Newton-Verfahren functions
# true for the normaler Newton-Verfahren and false for vereinfachte Newton-Verfahren
function setup_gauß(fh_F::Function,fh_DF::Function,a::Vector{Float64},bol::Bool)

          bol ? DF = fh_DF(a) : DF=V
          f  = -( fh_F(a))
          m = 0.0;
          al,ac = size(DF) # al :Ligne size  , ac : colomn size
          s = zeros(al)

          #Gauß Elimination
          for k= 1:(al-1)
              for i = (k+1):al
                  m = DF[i,k]/(DF[k,k])
                  DF[i,k] = 0
                  for j=(k+1):al
                      DF[i,j] = DF[i,j] - m*DF[k,j]
                  end
                    f[i]= f[i] - m*f[k]
              end
          end
          #Calculete Sk
          s[al] = f[al]/(DF[al,al])
          for k = (al-1):-1:1
              t = 0;
              for j = (k+1):al
                  t = t+DF[k,j]*s[j]
              end
              s[k]=(f[k]-t)/DF[k,k]
          end

          return s
end

# normaler newton verfahren
function Newton_Verfahren(x::Vector{Float64},bol::Bool)
  xneu = zeros(n)
  t=0 #Iteration
  while t < 1000

        xneu  = x + setup_gauß(model_F,model_DF,x,bol)

        temp = norm(xneu - x)
        if temp < epsilon_x && norm(model_F(xneu)) < epsilon_f
            bol ? print("Normaler Newton Verfahren succeed in :",t," Iterationen") : print("das vereinfachte Newton Verfahren succeed in :",t," Iterationen")

            return xneu
        end
        x = xneu
        t += 1
 end
  return false
end


#main
#initial
n = 60
a = 0
b = 1
h = (b-a)/n
#Abbruchkriterien
#x0 = fill(3.0,n)
x0 = fill(1.0,n)    # Startwert
#epsilon_x = 1e-17
epsilon_x = 1e-14   # Fehler-Kriterien
#epsilon_f = 1e-17
epsilon_f = 1e-13   # Residuums-Kriterium

V=model_DF(x0)       #we Fix DF for the Vereinfachte Newton Verfahren function in x0

u1 = Newton_Verfahren(x0,true)
u2 = Newton_Verfahren(x0,false)


# plot solution
tx = collect(0:1/n:1)


import Winston # include plot package

pl = Winston.plot(tx,u1,"b.")

Winston.setattr(pl,"title","IntegralGleichung")
Winston.setattr(pl,"xlabel","x")
Winston.setattr(pl,"ylabel","y")

Winston.display(pl)
