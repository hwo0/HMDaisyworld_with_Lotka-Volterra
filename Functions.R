## Supporting functions for daisyworld
# 7th Jun 2022
# Adjusted with Lotka-Volterra Equations (addition of Herbivores)
# Developed based on Dr.Robin Hankin's initial scripts
# Found at: https://github.com/RobinHankin/daisyworld

library("deSolve")

## First define physical parameters, all in one place:
parameters <- list(
    L  = 1,  # solar luminosity
    gamma = 0, # death rate
    parabolic = 0.003265,
    T_opt = 22.5,
    S     = 917, # SI!  In fig 1, W&L use ergs/cm^2s FFS
    sigma = 5.670374419e-8, # stefan's constant: SI
    triple_point = 273, # actually 273.15
    Ag = 0.50, # albedo of ground
    Aw = 0.75, # albedo of white daisies
    Ab = 0.25, # albedo of black daisies
    q = NA, 
    qdash = 20,    # caption in figure 1
    p = 1.0,  # proportion of fertile ground, fig 1
    K= 1, #rate herbivores eat daises
    C=0.8, #rate herbivores reproduce by eating daisies
    D=0.3  #rate herbivores die by natural cause 
)

## Next define R functions for the equations in Watson and Lovelock:

`bare_fertile` <- function(W,B){with(parameters,p-W-B)} #eq 3
`growth_rate` <- function(temp,parameters){with(parameters, 1-parabolic*(T_opt-temp)^2)} #eq 6
`T_eff` <- function(A,L,parameters){with(parameters, (S*L*(1-A)/sigma)^(0.25) - triple_point)} #eq 7
`albedo` <- function(W,B,parameters){with(parameters, (1-W-B)*Ag + W*Aw + B*Ab)} #eq 8
`T_white` <- function(T_e, A){with(parameters, qdash*(A-Aw) + T_e ) } #eq 9 (white)
`T_black` <- function(T_e, A){with(parameters, qdash*(A-Ab) + T_e ) } #eq 9 (black)


## Now define a function, watson(), for use with ode() 

`watson` <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
        x <- bare_fertile(W,B)
        A <- albedo(W,B,parameters)
        T_e <- T_eff(A,L,parameters)
        
        T_W <- T_white(T_e,A)
        T_B <- T_black(T_e,A)
        
        beta_W <- growth_rate(T_W,parameters)
        beta_B <- growth_rate(T_B,parameters)
        
        dW <-  W*(x*beta_W-gamma)-K*W*P
        dB <-  B*(x*beta_B-gamma)-K*B*P
        
        dP <- C*(W+B)*P - D*P 
        
        list(c(dW,dB,dP),T_e=T_e)
    })
}

