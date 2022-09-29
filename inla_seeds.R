data(Seeds)

formula <- r~x1 + x2 + x1*x2 + f(plate, model="iid")

model.regression <- inla(formula, data=Seeds,
                         family="binomial", 
                         Ntrials=n)
