import math
import numpy as np

### GEAR SET 1 ###

# INPUTS

power = 20  # horsepower
lifetime = 12000  # gear lifetime [hours]
inputSpeed = 1750  # rpm
outputSpeed = 85
minNumTeeth = 17  # From table xxxx
Ko = 2.00  # Adjustment factor from 9-1 -> Moderate Shock Input, Moderate Shock Output
serviceFactor = 1.00

# SOLVER

targetVelocityRatio = inputSpeed / outputSpeed
print("Target Gear Ratio: " + str(targetVelocityRatio))

print("GEAR SET 1 PARAMS \n\n\n")

targetVR = math.sqrt(targetVelocityRatio)  # Velocity Ratio for both sets of gears

pinionTeeth = minNumTeeth
gearTeeth = round(pinionTeeth * targetVR)
print("Gear Teeth: " + str(gearTeeth))

VR = gearTeeth / pinionTeeth
print("Velocity Ratio: " + str(VR))

adjustedPower = Ko * power
pitchDiameter = 7.00  # Pd: From 9-11, adjusted power corresponds to a pitch diameter of 7 in

gearDiameter = gearTeeth / pitchDiameter
pinionDiameter = pinionTeeth / pitchDiameter

adjustmentDiameter = 1 / pitchDiameter


print("Adjustment Diameter: " + str(adjustmentDiameter))
print("Gear Diameter: " + str(gearDiameter))
print("Pinion Diameter: " + str(pinionDiameter))

geartrainHeight = (3 / 2) * gearDiameter + pinionDiameter / 2 + 3 * adjustmentDiameter  # calculated Y (in)

print("Geartrain Height: " + str(geartrainHeight))

C = (pinionDiameter + gearDiameter) / (2 * pitchDiameter)  # center distance
tangentialVelocity = (math.pi * pitchDiameter * inputSpeed) / 12
transmittedLoad = 33000 * power / tangentialVelocity  # Wt

print("C: " + str(C))
print("Tangential velocity: " + str(tangentialVelocity))
print("wt: " + str(transmittedLoad))

nominalFaceWidth = 12 / pitchDiameter

print("Nominal Face Width: " + str(nominalFaceWidth))

# Material is steel
elasticCoefficient = 2300  # Cp
agmaQualityNumber = 6  # A6 - automotive transmission equiv. from 9-5

dynamicFactorsAt3207VT = {}  # Mapping between AGMA Quality Numbers and Dynamic Factors at Vt of 3207
dynamicFactorsAt3207VT[6] = 1.13
dynamicFactorsAt3207VT[7] = 1.23
dynamicFactorsAt3207VT[8] = 1.34
dynamicFactorsAt3207VT[9] = 1.46
dynamicFactorsAt3207VT[10] = 1.58
dynamicFactorsAt3207VT[11] = 1.75

dynamicFactor = dynamicFactorsAt3207VT[agmaQualityNumber]  # Kv

print("Dynamic Factor: " + str(dynamicFactor))

KMRatio = nominalFaceWidth / pinionDiameter

print("KM Ratio: " + str(KMRatio))

# Use this formula because 1 < F < 15 (See Figure 9-12)
pinionProportionFactor = nominalFaceWidth / (10 * pinionDiameter) - 0.0375 + 0.0125 * nominalFaceWidth  # from fig 9-12

print("Cpf: " + str(pinionProportionFactor))

# Assume commercial gears
meshAlignmentFactor = 0.127 + 0.0158 * nominalFaceWidth - (1.093 * pow(10, -4) * pow(nominalFaceWidth, 2))  # Cma

Km = 1 + meshAlignmentFactor + pinionProportionFactor

print("Mesh Alignment Factor: " + str(meshAlignmentFactor))

sizeFactor = 1  # Ks since pitch diameter is greater than 5 inches from Table 9-2

rimThicknessFactor = 1.0  # Kb since VR is greater than 1.2, thickness factor is 1.0

reliabilityFactor = 1.00  # Km From table 9-11, at 99% reliability

numPinionCycles = 60 * lifetime * inputSpeed  # Ncp
numGearCycles = 60 * lifetime * (inputSpeed / VR)  # Ncg

print("Number Pinion Cycles: " + str(numPinionCycles))
print("Number Gear Cycles: " + str(numGearCycles))


def bending_strength_stress_cycle_factor(numCycles):
    # This is linear because our load cycles is greater than 10^7 (Fig 9-21)
    return 1.3558 * numCycles ** -0.0178


def pitting_strength_stress_cycle_factor(numCycles):
    # This is linear because our load cycles is greater than 10^7 (Fig 9-22)
    return 1.4488 * numCycles ** -0.023


bendingStressCycleFactorPinion = bending_strength_stress_cycle_factor(numPinionCycles)  # Ynp from Figure 9-21
bendingStressCycleFactorGear = bending_strength_stress_cycle_factor(numGearCycles)  # Yng from Figure 9-21

pittingResistanceStressCycleFactorPinion = pitting_strength_stress_cycle_factor(numPinionCycles)  # Znp from Figure 9-22
pittingResistanceStressCycleFactorGear = pitting_strength_stress_cycle_factor(numGearCycles)  # Zng from Figure 9-22

print("Ynp: " + str(bendingStressCycleFactorPinion))
print("Yng: " + str(bendingStressCycleFactorGear))

print("Znp: " + str(pittingResistanceStressCycleFactorPinion))
print("Zng: " + str(pittingResistanceStressCycleFactorGear))

# From Figure 9-10: Np = 17, Ng = 77
bendingGeometryFactorPinion = 0.295  # Jp from Fig 9-17
bendingGeometryFactorGear = 0.41  # Jg from Fig 9-17

# From Figure 9-17 VR = 4.54 Np = 17

pittingGeometryFactor = 0.105  # I from Figure 9-17, using 20deg pressure angle

allowableBendingStressNumberPinion = (((transmittedLoad * pitchDiameter) / (nominalFaceWidth * bendingGeometryFactorPinion)) * Ko \
                                      * Km * sizeFactor * rimThicknessFactor * dynamicFactor) * (
                                                 serviceFactor * reliabilityFactor / (
                                             bendingStressCycleFactorPinion))  # Satp
allowableBendingStressNumberGear = (((transmittedLoad * pitchDiameter) / (nominalFaceWidth * bendingGeometryFactorGear)) * Ko \
                                    * Km * sizeFactor * rimThicknessFactor * dynamicFactor) * (
                                               serviceFactor * reliabilityFactor / (
                                           bendingStressCycleFactorGear))  # Satg

allowableContactStressNumberPinion = ((elasticCoefficient * reliabilityFactor * serviceFactor) /
                                      pittingResistanceStressCycleFactorPinion) * ((transmittedLoad * Ko * sizeFactor * Km *
                                                                                    dynamicFactor) / (nominalFaceWidth *
                                                                                                      pinionDiameter *
                                                                                                      pittingGeometryFactor)) ** 0.5  # Sacp

allowableContactStressNumberGear = ((elasticCoefficient * reliabilityFactor * serviceFactor) /
                                      pittingResistanceStressCycleFactorGear) * ((transmittedLoad * Ko * sizeFactor * Km *
                                                                                  dynamicFactor) / (nominalFaceWidth *
                                                                                                      gearDiameter *
                                                                                                      pittingGeometryFactor)) ** 0.5  # Sacg

print("Bending Stress Number Pinion: " + str(allowableBendingStressNumberPinion))
print("Bending Stress Number Gear: " + str(allowableBendingStressNumberGear))
print("Contact Stress Number Pinion: " + str(allowableContactStressNumberPinion))
print("Contact Stress Number Gear: " + str(allowableContactStressNumberGear))

limitingStress = max([allowableBendingStressNumberPinion, allowableBendingStressNumberGear, allowableContactStressNumberGear, allowableContactStressNumberPinion])

print("Limiting Factor: " + str(limitingStress))

def calc_brinell_hardness(limitingStress):
    return (limitingStress / 1000 - 29.10) / 0.322

brinellHardness = calc_brinell_hardness(limitingStress) # HB from Fig 9-19 using Grade 1 Steel
print("Brinell Hardness: " + str(brinellHardness))

# Use appendix 3 to select a steel -> Use 1144 Cold-Drawn (HB = 200)


### GEAR SET 2 ###

print("\n\n\n GEAR SET 2 PARAMS \n\n\n ")

targetVR = math.sqrt(targetVelocityRatio)  # Velocity Ratio for both sets of gears

pinionTeeth = minNumTeeth
gearTeeth = round(pinionTeeth * targetVR)
print("Gear Teeth: " + str(gearTeeth))

VR = gearTeeth / pinionTeeth
print("Velocity Ratio: " + str(VR))

inputSpeed = inputSpeed / VR

adjustedPower = Ko * power
pitchDiameter = 7.00  # Pd: From 9-11, adjusted power corresponds to a pitch diameter of 7 in

gearDiameter = gearTeeth / pitchDiameter
pinionDiameter = pinionTeeth / pitchDiameter

adjustmentDiameter = 1 / pitchDiameter

print("Gear Diameter: " + str(gearDiameter))
print("Pinion Diameter: " + str(pinionDiameter))

geartrainHeight = (3 / 2) * gearDiameter + pinionDiameter / 2 + 3 * adjustmentDiameter  # calculated Y (in)

print("Geartrain Height: " + str(geartrainHeight))

C = (pinionDiameter + gearDiameter) / (2 * pitchDiameter)  # center distance
tangentialVelocity = (math.pi * pitchDiameter * inputSpeed) / 12
transmittedLoad = 33000 * power / tangentialVelocity  # Wt

print("C: " + str(C))
print("Tangential velocity: " + str(tangentialVelocity))
print("wt: " + str(transmittedLoad))

nominalFaceWidth = 12 / pitchDiameter

print("Nominal Face Width: " + str(nominalFaceWidth))

# Material is steel
elasticCoefficient = 2300  # Cp
agmaQualityNumber = 6  # A6 - automotive transmission equiv. from 9-5

dynamicFactorsAt3207VT = {}  # Mapping between AGMA Quality Numbers and Dynamic Factors at Vt of 708
dynamicFactorsAt3207VT[6] = 1.07
dynamicFactorsAt3207VT[7] = 1.12
dynamicFactorsAt3207VT[8] = 1.18
dynamicFactorsAt3207VT[9] = 1.24
dynamicFactorsAt3207VT[10] = 1.31
dynamicFactorsAt3207VT[11] = 1.39

dynamicFactor = dynamicFactorsAt3207VT[agmaQualityNumber]  # Kv

print("Dynamic Factor: " + str(dynamicFactor))

KMRatio = nominalFaceWidth / pinionDiameter

print("KM Ratio: " + str(KMRatio))

# Use this formula because 1 < F < 15 (See Figure 9-12)
pinionProportionFactor = nominalFaceWidth / (10 * pinionDiameter) - 0.0375 + 0.0125 * nominalFaceWidth  # from fig 9-12

print("Cpf: " + str(pinionProportionFactor))

# Assume commercial gears
meshAlignmentFactor = 0.127 + 0.0158 * nominalFaceWidth - (1.093 * pow(10, -4) * pow(nominalFaceWidth, 2))  # Cma

Km = 1 + meshAlignmentFactor + pinionProportionFactor

print("Mesh Alignment Factor: " + str(meshAlignmentFactor))

sizeFactor = 1  # Ks since pitch diameter is greater than 5 inches from Table 9-2

rimThicknessFactor = 1.0  # Kb since VR is greater than 1.2, thickness factor is 1.0

reliabilityFactor = 1.00  # Km From table 9-11, at 99% reliability

numPinionCycles = 60 * lifetime * inputSpeed  # Ncp
numGearCycles = 60 * lifetime * (inputSpeed / VR)  # Ncg

print("Number Pinion Cycles: " + str(numPinionCycles))
print("Number Gear Cycles: " + str(numGearCycles))

bendingStressCycleFactorPinion = bending_strength_stress_cycle_factor(numPinionCycles)  # Ynp from Figure 9-21
bendingStressCycleFactorGear = bending_strength_stress_cycle_factor(numGearCycles)  # Yng from Figure 9-21

pittingResistanceStressCycleFactorPinion = pitting_strength_stress_cycle_factor(numPinionCycles)  # Znp from Figure 9-22
pittingResistanceStressCycleFactorGear = pitting_strength_stress_cycle_factor(numGearCycles)  # Zng from Figure 9-22

print("Ynp: " + str(bendingStressCycleFactorPinion))
print("Yng: " + str(bendingStressCycleFactorGear))

print("Znp: " + str(pittingResistanceStressCycleFactorPinion))
print("Zng: " + str(pittingResistanceStressCycleFactorGear))

# From Figure 9-10: Np = 17, Ng = 77
bendingGeometryFactorPinion = 0.295  # Jp from Fig 9-17
bendingGeometryFactorGear = 0.41  # Jg from Fig 9-17

# From Figure 9-17 VR = 4.54 Np = 17

pittingGeometryFactor = 0.105  # I from Figure 9-17, using 20deg pressure angle

allowableBendingStressNumberPinion = (((transmittedLoad * pitchDiameter) / (nominalFaceWidth * bendingGeometryFactorPinion)) * Ko \
                                      * Km * sizeFactor * rimThicknessFactor * dynamicFactor) * (
                                                 serviceFactor * reliabilityFactor / (
                                             bendingStressCycleFactorPinion))  # Satp
allowableBendingStressNumberGear = (((transmittedLoad * pitchDiameter) / (nominalFaceWidth * bendingGeometryFactorGear)) * Ko \
                                    * Km * sizeFactor * rimThicknessFactor * dynamicFactor) * (
                                               serviceFactor * reliabilityFactor / (
                                           bendingStressCycleFactorGear))  # Satg

allowableContactStressNumberPinion = ((elasticCoefficient * reliabilityFactor * serviceFactor) /
                                      pittingResistanceStressCycleFactorPinion) * ((transmittedLoad * Ko * sizeFactor * Km *
                                                                                    dynamicFactor) / (nominalFaceWidth *
                                                                                                      pinionDiameter *
                                                                                                      pittingGeometryFactor)) ** 0.5  # Sacp

allowableContactStressNumberGear = ((elasticCoefficient * reliabilityFactor * serviceFactor) /
                                      pittingResistanceStressCycleFactorGear) * ((transmittedLoad * Ko * sizeFactor * Km *
                                                                                  dynamicFactor) / (nominalFaceWidth *
                                                                                                      gearDiameter *
                                                                                                      pittingGeometryFactor)) ** 0.5  # Sacg

print("Bending Stress Number Pinion: " + str(allowableBendingStressNumberPinion))
print("Bending Stress Number Gear: " + str(allowableBendingStressNumberGear))
print("Contact Stress Number Pinion: " + str(allowableContactStressNumberPinion))
print("Contact Stress Number Gear: " + str(allowableContactStressNumberGear))

limitingStress = max([allowableBendingStressNumberPinion, allowableBendingStressNumberGear, allowableContactStressNumberGear, allowableContactStressNumberPinion])

print("Limiting Factor: " + str(limitingStress))

brinellHardness = calc_brinell_hardness(limitingStress) # HB from Fig 9-19 using Grade 1 Steel
print("Brinell Hardness: " + str(brinellHardness))

print("Output Speed: " + str(inputSpeed / VR))
# Use appendix 3 to select a steel -> Use SAE 4150 OQT 700 (HB = 495)

