#using Plots
using DataFrames
using Gadfly

function SIR2pop(βₒ, βᵤ; baseMortalityRatio::Float64=1/30, ageMortalityRatio::Float64 = 1/50, N::Float64 = 300e6, sf::Float64 = .25, stepUp::Int64 = 180, overCrowdCost::Int64 = 2, overCrowdThresh::Float64 = .5e6, nStep::Int64 = 1000, stepSize::Int64 = 1)

  #All variables and parameters come with an "o" and "y" subscript to indicate the "Over 65" and "under 65" ("Young") populations.

  #hard coded parameters"
  over = .15 #fraction of population over-65,
  Iₒ = 1000 #number of infected individuals
  Iᵤ= 1000
  year = 365 #year, time point at which transmission levels return to normal values
  revertR = 2.8 #COVID-19 R0 (the average number of individuals infected by a single individual in a completely susceptible population)"
  recoveryTime = 14 #recovery time in days. Taken to be 14 day
  hmratio = 10 #the hospitalization/mortality ratio

  #population variables
  Nₒ0 = N  *  over
  Nᵤ0 = N - Nₒ0

  #setting up initial conditions
  Sₒ = Nₒ0  *  sf #number susceptible
  Sᵤ = Nᵤ0  *  sf
  Rₒ = Nₒ0 - Sₒ #number recovered
  Rᵤ = Nᵤ0 - Sᵤ
  Mᵤ = Hᵤ = 0 #number of mortalities and hospitalizations is 0 at the beginning
  Mₒ = Hₒ = 0

  #default values for gamma and delta variables that specify the rates of leaving the infected state
  gₒ = gᵤ = 1 / recoveryTime #recovery rate
  Dₒ = gₒ * baseMortalityRatio #mortality rate for the "over 65" population
  Dᵤ = (gᵤ * baseMortalityRatio) * ageMortalityRatio #mortality rate for the "under 65" population,

  #define alpha, the total rate of leaving the infected state
  αₒ = gₒ + Dₒ  #Recovery rate  +  Death rate for old people
  αᵤ = gᵤ + Dᵤ  #Recovery rate  +  Death rate for young people

  #set the actual beta parameters
  βₒₒ = βₒ * (αₒ * Nₒ0 / N + αᵤ * Nᵤ0 / N) / sf
  βᵤᵤ = βᵤ * (αₒ * Nₒ0 / N + αᵤ * Nᵤ0 / N) / sf

  #Return minimum between boo and byy
  βᵤₒ = min(βₒₒ, βᵤᵤ) #we assume that the cross age group transmission is the minimum of the two

  out = zeros(nStep + 1, 10)
  colnames = ["Sₒ" "Iₒ" "Rₒ" "Hₒ" "Mₒ" "Sᵤ" "Iᵤ " "Rᵤ" "Hᵤ" "Mᵤ"]

  #compute the final betas's. These are the betas that correspond to natural transmission levels without intervention that the model reverts to at 365 days
  βₒₒFinal = revertR * (αₒ * Nₒ0 / N + αᵤ * Nᵤ0 / N) / sf

  βᵤᵤFinal = revertR * (αₒ * Nₒ0 / N + αᵤ * Nᵤ0 / N) / sf
  βᵤₒFinal = min(βₒₒFinal,βᵤᵤFinal)

  #save the difference between final and initial values so we can step up
  βₒₒDiff = βₒₒFinal - βₒₒ
  βᵤᵤDiff = βᵤᵤFinal - βᵤᵤ
  βᵤₒDiff = βᵤₒFinal - βᵤₒ

  for s in 1:nStep
    #println(s)
    totalHospital = Hₒ + Hᵤ #Number of hospitalizations for old and young
    if totalHospital > overCrowdThresh
      overFrac = (totalHospital - overCrowdThresh) / totalHospital #the ratio of hospitalized individuals that are over capacity
      multiplier = overFrac * overCrowdCost + (1 - overFrac) * 1 # a weighted average of mortality for the below capacity and above
      #capacity fraction of hospitalized individuals
      thisDₒ  =  multiplier * Dₒ
      thisDᵤ = multiplier * Dᵤ
    else
      thisDₒ = Dₒ
      thisDᵤ = Dᵤ
    end
    thisαₒ = gₒ  +  thisDₒ
    thisαᵤ = gᵤ  +  thisDᵤ

    out[s,:] = [Sₒ Iₒ Rₒ Hₒ Mₒ Sᵤ Iᵤ  Rᵤ Hᵤ Mᵤ]

    newIᵤ  = Iᵤ  + (βᵤᵤ * Sᵤ * Iᵤ  / N + βᵤₒ * Sᵤ * Iₒ / N - thisαᵤ * Iᵤ ) * stepSize
    Sᵤ = Sᵤ - (βᵤᵤ * Sᵤ * Iᵤ  / N + βᵤₒ * Sᵤ * Iₒ / N) * stepSize
    Rᵤ = Rᵤ + gᵤ * Iᵤ  * stepSize
    Mᵤ = Mᵤ + thisDᵤ * Iᵤ  * stepSize
    Hᵤ = hmratio * recoveryTime * Dᵤ * Iᵤ  #dy instead of thisdy is used so that the excess mortalities (due to overcrowding) doesn't creat extra hospitalizations

    newIₒ = Iₒ + (βᵤₒ * Sₒ * Iᵤ  / N + βₒₒ * Sₒ * Iₒ / N - thisαₒ * Iₒ) * stepSize
    Sₒ = Sₒ - (βᵤₒ * Sₒ * Iᵤ  / N + βₒₒ * Sₒ * Iₒ / N) * stepSize
    Rₒ = Rₒ + gₒ * Iₒ * stepSize
    Mₒ = Mₒ + thisDₒ * Iₒ * stepSize
    Hₒ = hmratio * recoveryTime * Dₒ * Iₒ

    Iₒ = newIₒ
    Iᵤ  = newIᵤ

    #assumption number 2: mitigation strategies must survive reintroduction
    Iₒ = max(Iₒ,1)
    Iᵤ  = max(Iᵤ ,1)

    #step up infection rates to normal levels by one year starting at step.up
    if s * stepSize > stepUp && s * stepSize <= year
      βₒₒ = βₒₒ + βₒₒDiff / ((year - stepUp) / stepSize)
      βᵤᵤ = βᵤᵤ + βᵤᵤDiff / ((year - stepUp) / stepSize)
      βᵤₒ = βᵤₒ + βᵤₒDiff / ((year - stepUp) / stepSize)
    end
  #End of for loop
  end
  out[nStep  +  1, :] = [Sₒ Iₒ Rₒ Hₒ Mₒ Sᵤ Iᵤ Rᵤ Hᵤ Mᵤ]
  println([Mᵤ Mₒ Mᵤ + Mₒ])

  return out
end

#A is the result of STR2pop function
A = SIR2pop(1, 1)
#Converts array matrix to DataFrame
B = convert(DataFrame, A)
#Adds column names
rename!(B,Dict(:x1 => :Sₒ, :x2 => :Iₒ, :x3 => :Rₒ, :x4 => :Hₒ, :x5 => :Mₒ, :x6 => :Sᵤ, :x7 => :Iᵤ , :x8 => :Rᵤ, :x9 => :Hᵤ, :x10 => :Mᵤ))
#Display B to REPL
show(B)
