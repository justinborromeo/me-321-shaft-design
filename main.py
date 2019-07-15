import math

# INPUTS

power = 20 # horsepower
lifetime = 12000 # gear lifetime [hours]
inputSpeed = 1750 # rpm
outputSpeed = 85
minNumTeeth = 17 # From table xxxx
Ko = 2.00 # Adjustment factor from 9-1 -> Moderate Shock Input, Moderate Shock Output


# SOLVER

targetVelocityRatio = inputSpeed / outputSpeed
print("Target Gear Ratio: " + str(targetVelocityRatio))

VR = math.sqrt(targetVelocityRatio) # Velocity Ratio for both sets of gears
print("Velocity Ratio: " + str(VR))

pinionTeeth = minNumTeeth
gearTeeth = round(pinionTeeth * VR)
print("Gear Teeth: " + str(gearTeeth))

velocityRatio = gearTeeth/pinionTeeth

adjustedPower = Ko * power
pitchDiameter = 7.00 # Pd: From 9-11, adjusted power corresponds to a pitch diameter of 7 in

gearDiameter = gearTeeth / pitchDiameter
pinionDiameter = pinionTeeth / pitchDiameter

adjustmentDiameter = 1 / pitchDiameter

print("Gear Diameter: " + str(gearDiameter))
print("Pinion Diameter: " + str(pinionDiameter))

geartrainHeight = (3/2) * gearDiameter + pinionDiameter / 2 + 3 * adjustmentDiameter # calculated Y (in)

print("Geartrain Height: " + str(geartrainHeight))

C = (pinionDiameter + gearDiameter) / (2 * pitchDiameter) # center distance
tangentialVelocity = (math.pi * pitchDiameter * inputSpeed) / 12
wt = 33000 * power / tangentialVelocity # TODO what is this?

print("C: " + str(C))
print("Tangential velocity: " + str(tangentialVelocity))
print("wt: " + str(wt))

nominalFaceWidth = 12 / pitchDiameter

print("Nominal Face Width: " + str(nominalFaceWidth))

# Material is steel
elasticCoefficient = 2300 # Cp
agmaQualityNumber = 6 # A6 - automotive transmission equiv. from 9-5

dynamicFactorsAt3207VT = {} # Mapping between AGMA Quality Numbers and Dynamic Factors at Vt of 3207
dynamicFactorsAt3207VT[6] = 1.13
dynamicFactorsAt3207VT[7] = 1.23
dynamicFactorsAt3207VT[8] = 1.34
dynamicFactorsAt3207VT[9] = 1.46
dynamicFactorsAt3207VT[10] = 1.58
dynamicFactorsAt3207VT[11] = 1.75

dynamicFactor = dynamicFactorsAt3207VT[agmaQualityNumber] # Kv

print("Dynamic Factor: " + str(dynamicFactor))

KMRatio = nominalFaceWidth / pinionDiameter

print("KM Ratio: " + str(KMRatio))

# Use this formula because 1 < F < 15 (See Figure 9-12)
pinionProportionFactor = nominalFaceWidth / (10 * pinionDiameter) - 0.0375 + 0.0125 * nominalFaceWidth # from fig 9-12

print("Cpf: " + str(pinionProportionFactor))

# Assume commercial gears
meshAlignmentFactor = 0.127 + 0.0158 * nominalFaceWidth - (1.093 * pow(10, -4) * pow(nominalFaceWidth, 2)) # Cma

print("Mesh Alignment Factor" + str(meshAlignmentFactor))

sizeFactor = 1 # since pitch diameter is greater than 5 inches from Table 9-2

rimThicknessFactor = 1.0 # since VR is greater than 1.2, thickness factor is 1.0

reliabilityFactor = 1.00 # From table 9-11, at 99% reliability

numPinionCycles = 60 * lifetime * inputSpeed # Ncp
numGearCycles = 60 * lifetime * (inputSpeed / velocityRatio) # Ncg

print("Number Pinion Cycles: " + str(numPinionCycles))
print("Number Gear Cycles: " + str(numGearCycles))




