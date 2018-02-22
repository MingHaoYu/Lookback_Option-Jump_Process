# Lookback_Option-Jump_Process
Lookback option: <br>
Pricing Fixed Strike Lookback Call and Put Options.

parameters: <br>
r=0.03, stock price=98, strike price=100, volatility= 0.12 to 0.48 <br>
Call option payoff function = max(Smax-k, 0) <br>
Put option payoff function = max(k-Smin, 0) <br>


Jump-diffusion process: <br>
dVt = Vt * mu * dt + Vt * sigma * dWt + Vt * gamma * dJt <br>
where mu, sigma, gamma <0, and V0 are given, J is a Poisson process, with intensity lambda1, independent of the Brownian Motion process W.
