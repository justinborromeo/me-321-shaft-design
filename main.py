import math
import numpy as np
import matplotlib.pyplot as plt

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

print("GEAR SET 1 PARAMS \n")
print("1) Target Velocity Ratio: " + str(targetVelocityRatio))
targetVR = math.sqrt(targetVelocityRatio)  # Velocity Ratio for both sets of gears
print("2) Target Gearset Velocity Ratio:" + str(targetVR))
print("3) Num Pinion Teeth: " + str(minNumTeeth))
pinionTeeth = minNumTeeth
gearTeeth = round(pinionTeeth * targetVR)
print("4) Num Gear Teeth: " + str(gearTeeth))

VR = gearTeeth / pinionTeeth
print("5) Actual Velocity Ratio: " + str(VR))

adjustedPower = Ko * power
print("6) Design Power [W]: " + str(adjustedPower))
pitchDiameter = 7.00  # Pd: From 9-11, adjusted power corresponds to a pitch diameter of 7 in
print("7) Pitch Diameter [in]: " + str(pitchDiameter))
gearDiameter = gearTeeth / pitchDiameter
pinionDiameter = pinionTeeth / pitchDiameter

adjustmentDiameter = 1 / pitchDiameter

print("8) Gear Diameter [in]: " + str(gearDiameter))
print("8) Pinion Diameter [in]: " + str(pinionDiameter))

addendum = adjustmentDiameter
dedendum = 1.25 / pitchDiameter

print("Addendum: " + str(addendum))
print("Dedendum: " + str(dedendum))

addendumDiameter = 2 * addendum + pinionDiameter
dedendumDiameter = pinionDiameter - 2 * dedendum

print("Pinion Addendum Diameter: " + str(addendumDiameter))
print("Pinion Dedendum Diameter: " + str(dedendumDiameter))

addendumDiameter = 2 * addendum + gearDiameter
dedendumDiameter = gearDiameter - 2 * dedendum

print("Gear Addendum Diameter: " + str(addendumDiameter))
print("Gear Dedendum Diameter: " + str(dedendumDiameter))

geartrainHeight = (3 / 2) * gearDiameter + pinionDiameter / 2 + 3 * adjustmentDiameter  # calculated Y (in)

print("9) Geartrain Height [in]: " + str(geartrainHeight))

C = (pinionDiameter + gearDiameter) / (2 * pitchDiameter)  # center distance
print("10) Center Distance [in]: " + str(C))
tangentialVelocity = (math.pi * pitchDiameter * inputSpeed) / 12
transmittedLoad = 33000 * power / tangentialVelocity  # Wt

print("11) Pitch Line Speed [ft/min]: " + str(tangentialVelocity))
print("12) Transmitted Load: " + str(transmittedLoad))

nominalFaceWidth = 12 / pitchDiameter

print("13) Nominal Face Width: " + str(nominalFaceWidth))

# Material is steel
elasticCoefficient = 2300  # Cp
print("14) Elastic Coefficient: " + str(elasticCoefficient))

agmaQualityNumber = 6  # A6 - automotive transmission equiv. from 9-5
print("15) AGMA Quality Number: " + str(agmaQualityNumber))

dynamicFactorsAt3207VT = {}  # Mapping between AGMA Quality Numbers and Dynamic Factors at Vt of 3207
dynamicFactorsAt3207VT[6] = 1.13
dynamicFactorsAt3207VT[7] = 1.23
dynamicFactorsAt3207VT[8] = 1.34
dynamicFactorsAt3207VT[9] = 1.46
dynamicFactorsAt3207VT[10] = 1.58
dynamicFactorsAt3207VT[11] = 1.75

dynamicFactor = dynamicFactorsAt3207VT[agmaQualityNumber]  # Kv

print("16) Dynamic Factor: " + str(dynamicFactor))

KMRatio = nominalFaceWidth / pinionDiameter

print("17) F/Dp Ratio: " + str(KMRatio))

# Use this formula because 1 < F < 15 (See Figure 9-12)
pinionProportionFactor = nominalFaceWidth / (10 * pinionDiameter) - 0.0375 + 0.0125 * nominalFaceWidth  # from fig 9-12

print("18) Pinion Proportion Factor: " + str(pinionProportionFactor))

# Assume commercial gears
meshAlignmentFactor = 0.127 + 0.0158 * nominalFaceWidth - (1.093 * pow(10, -4) * pow(nominalFaceWidth, 2))  # Cma
print("19) Mesh Alignment Factor: " + str(meshAlignmentFactor))

Km = 1 + meshAlignmentFactor + pinionProportionFactor

print("20) Km: " + str(Km))

sizeFactor = 1  # Ks since pitch diameter is greater than 5 inches from Table 9-2
print("21) Ks: " + str(sizeFactor))

rimThicknessFactor = 1.0  # Kb since VR is greater than 1.2, thickness factor is 1.0
print("22) Rim Thickness Factor: " + str(rimThicknessFactor))

reliabilityFactor = 1.00  # Km From table 9-11, at 99% reliability
print("23) Reliability Factor: " + str(reliabilityFactor))

numPinionCycles = 60 * lifetime * inputSpeed  # Ncp
numGearCycles = 60 * lifetime * (inputSpeed / VR)  # Ncg

print("24) Number Pinion Cycles: " + str(numPinionCycles))
print("24) Number Gear Cycles: " + str(numGearCycles))


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

print("25) Ynp: " + str(bendingStressCycleFactorPinion))
print("25) Yng: " + str(bendingStressCycleFactorGear))

print("26) Znp: " + str(pittingResistanceStressCycleFactorPinion))
print("26) Zng: " + str(pittingResistanceStressCycleFactorGear))

# From Figure 9-10: Np = 17, Ng = 77
bendingGeometryFactorPinion = 0.295  # Jp from Fig 9-17
bendingGeometryFactorGear = 0.41  # Jg from Fig 9-17
print("27) Bending Geometry Factor Pinion " + str(bendingGeometryFactorPinion))
print("27) Bending Geometry Factor Gear " + str(bendingGeometryFactorGear))

# From Figure 9-17 VR = 4.54 Np = 17

pittingGeometryFactor = 0.105  # I from Figure 9-17, using 20deg pressure angle

print("28) Pitting Geometry Factor: " + str(pittingGeometryFactor))

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

print("29) Bending Stress Number Pinion: " + str(allowableBendingStressNumberPinion))
print("29) Bending Stress Number Gear: " + str(allowableBendingStressNumberGear))
print("30) Contact Stress Number Pinion: " + str(allowableContactStressNumberPinion))
print("30) Contact Stress Number Gear: " + str(allowableContactStressNumberGear))

limitingStress = max([allowableBendingStressNumberPinion, allowableBendingStressNumberGear, allowableContactStressNumberGear, allowableContactStressNumberPinion])

print("31) Limiting Stress: " + str(limitingStress))

def calc_brinell_hardness(limitingStress):
    return (limitingStress / 1000 - 29.10) / 0.322

brinellHardness = calc_brinell_hardness(limitingStress) # HB from Fig 9-19 using Grade 1 Steel
print("32) Brinell Hardness: " + str(brinellHardness))

# Use appendix 3 to select a steel -> Use 1144 Cold-Drawn (HB = 200)
print("33) Choose 1144 Cold-Drawn Steel with HB = 200")

### GEAR SET 2 ###

print("\nGEAR SET 2 PARAMS \n")

print("1) Target Velocity Ratio: " + str(targetVR))
targetVR = math.sqrt(targetVelocityRatio)  # Velocity Ratio for both sets of gears
print("2) Target Gearset Velocity Ratio:" + str(targetVR))
pinionTeeth = minNumTeeth
print("3) Num Pinion Teeth: " + str(minNumTeeth))
gearTeeth = round(pinionTeeth * targetVR)
print("4) Num Gear Teeth: " + str(gearTeeth))

VR = gearTeeth / pinionTeeth
print("5) Actual Velocity Ratio: " + str(VR))

inputSpeed = inputSpeed / VR

adjustedPower = Ko * power
print("6) Design Power: " + str(adjustedPower))
pitchDiameter = 7.00  # Pd: From 9-11, adjusted power corresponds to a pitch diameter of 7 in
print("7) Pitch Diameter: " + str(pitchDiameter))

gearDiameter = gearTeeth / pitchDiameter
pinionDiameter = pinionTeeth / pitchDiameter

adjustmentDiameter = 1 / pitchDiameter

print("8) Gear Diameter: " + str(gearDiameter))
print("8) Pinion Diameter: " + str(pinionDiameter))

addendumDiameter = 2 * addendum + pinionDiameter
dedendumDiameter = pinionDiameter - 2 * dedendum

print("Pinion Addendum Diameter: " + str(addendumDiameter))
print("Pinion Dedendum Diameter: " + str(dedendumDiameter))

addendumDiameter = 2 * addendum + gearDiameter
dedendumDiameter = gearDiameter - 2 * dedendum

print("Gear Addendum Diameter: " + str(addendumDiameter))
print("Gear Dedendum Diameter: " + str(dedendumDiameter))

geartrainHeight = (3 / 2) * gearDiameter + pinionDiameter / 2 + 3 * adjustmentDiameter  # calculated Y (in)

print("9) Geartrain Height: " + str(geartrainHeight))

C = (pinionDiameter + gearDiameter) / (2 * pitchDiameter)  # center distance
print("10) Center Distance: " + str(C))
tangentialVelocity = (math.pi * pitchDiameter * inputSpeed) / 12
transmittedLoad = 33000 * power / tangentialVelocity  # Wt

print("11) Pitch Line Speed: " + str(tangentialVelocity))
print("12) Transmitted Load: " + str(transmittedLoad))

nominalFaceWidth = 12 / pitchDiameter

print("13) Nominal Face Width: " + str(nominalFaceWidth))

# Material is steel
elasticCoefficient = 2300  # Cp
print("14) Elastic Coefficient: " + str(elasticCoefficient))
agmaQualityNumber = 6  # A6 - automotive transmission equiv. from 9-5
print("15) AGMA Quality Number: " + str(agmaQualityNumber))

# Mapping between AGMA Quality Numbers and Dynamic Factors at Vt of 708
dynamicFactorsAt3207VT = {6: 1.07, 7: 1.12, 8: 1.18, 9: 1.24, 10: 1.31, 11: 1.39}

dynamicFactor = dynamicFactorsAt3207VT[agmaQualityNumber]  # Kv

print("16) Dynamic Factor: " + str(dynamicFactor))

KMRatio = nominalFaceWidth / pinionDiameter

print("17) F/Dp Ratio: " + str(KMRatio))

# Use this formula because 1 < F < 15 (See Figure 9-12)
pinionProportionFactor = nominalFaceWidth / (10 * pinionDiameter) - 0.0375 + 0.0125 * nominalFaceWidth  # from fig 9-12

print("18) Pinion Proportion Factor: " + str(pinionProportionFactor))

# Assume commercial gears
meshAlignmentFactor = 0.127 + 0.0158 * nominalFaceWidth - (1.093 * pow(10, -4) * pow(nominalFaceWidth, 2))  # Cma

print("19) Mesh Alignment Factor: " + str(meshAlignmentFactor))

Km = 1 + meshAlignmentFactor + pinionProportionFactor

print("20) Km: " + str(Km))

sizeFactor = 1  # Ks since pitch diameter is greater than 5 inches from Table 9-2

print("21) Ks: " + str(sizeFactor))

rimThicknessFactor = 1.0  # Kb since VR is greater than 1.2, thickness factor is 1.0

print("22) Rim Thickness Factor: " + str(rimThicknessFactor))

reliabilityFactor = 1.00  # Km From table 9-11, at 99% reliability

print("23) Reliability Factor: " + str(reliabilityFactor))

numPinionCycles = 60 * lifetime * inputSpeed  # Ncp
numGearCycles = 60 * lifetime * (inputSpeed / VR)  # Ncg

print("24) Number Pinion Cycles: " + str(numPinionCycles))
print("24) Number Gear Cycles: " + str(numGearCycles))

bendingStressCycleFactorPinion = bending_strength_stress_cycle_factor(numPinionCycles)  # Ynp from Figure 9-21
bendingStressCycleFactorGear = bending_strength_stress_cycle_factor(numGearCycles)  # Yng from Figure 9-21

pittingResistanceStressCycleFactorPinion = pitting_strength_stress_cycle_factor(numPinionCycles)  # Znp from Figure 9-22
pittingResistanceStressCycleFactorGear = pitting_strength_stress_cycle_factor(numGearCycles)  # Zng from Figure 9-22

print("25) Ynp: " + str(bendingStressCycleFactorPinion))
print("25) Yng: " + str(bendingStressCycleFactorGear))

print("26) Znp: " + str(pittingResistanceStressCycleFactorPinion))
print("26) Zng: " + str(pittingResistanceStressCycleFactorGear))

# From Figure 9-10: Np = 17, Ng = 77
bendingGeometryFactorPinion = 0.295  # Jp from Fig 9-17
bendingGeometryFactorGear = 0.41  # Jg from Fig 9-17
print("27) Bending Geometry Factor Pinion " + str(bendingGeometryFactorPinion))
print("27) Bending Geometry Factor Gear " + str(bendingGeometryFactorGear))

# From Figure 9-17 VR = 4.54 Np = 17

pittingGeometryFactor = 0.105  # I from Figure 9-17, using 20deg pressure angle
print("28) Pitting Geometry Factor: " + str(pittingGeometryFactor))
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

print("29) Bending Stress Number Pinion: " + str(allowableBendingStressNumberPinion))
print("29) Bending Stress Number Gear: " + str(allowableBendingStressNumberGear))
print("30) Contact Stress Number Pinion: " + str(allowableContactStressNumberPinion))
print("30) Contact Stress Number Gear: " + str(allowableContactStressNumberGear))

limitingStress = max([allowableBendingStressNumberPinion, allowableBendingStressNumberGear, allowableContactStressNumberGear, allowableContactStressNumberPinion])

print("31) Limiting Stress: " + str(limitingStress))

brinellHardness = calc_brinell_hardness(limitingStress) # HB from Fig 9-19 using Grade 1 Steel
print("32) Brinell Hardness: " + str(brinellHardness))

# Use appendix 3 to select a steel -> Use SAE 4150 OQT 700 (HB = 495)
print("33) Choose SAE 4150 OQT 700 with HB = 495")

print("Final Output Speed: " + str(inputSpeed / VR))


# Shaft 1

print("\nSHAFT 1\n")

powerIn = 20 * 63000 # watts
lengthSec1 = 0.5 # inches
lengthSec2 = 1
lengthSec3ToMiddle = 1.107
lengthSec3MiddleToEnd = 0.857
lengthSec4 = 1
lengthSec5 = 0.25
lengthSec6 = 0.25

lengthAToCenter = lengthSec1 + lengthSec2 + lengthSec3ToMiddle
lengthCenterToC = lengthSec3MiddleToEnd + lengthSec4 + lengthSec5 + lengthSec6;
totalLength = lengthAToCenter + lengthCenterToC

def plot_shear_diagram(force, lengthAToCenter, lengthCenterToC, figureNumber, plotTitle):
    rc = force * lengthAToCenter / (lengthAToCenter + lengthCenterToC)

    ra = force - rc
    print("10) Reaction Force C [lbs]: " + str(rc))
    print("10) Reaction Force A [lbs]: " + str(ra))
    x = np.linspace(0, lengthAToCenter + lengthCenterToC, discreteSubdivisions)
    y = np.linspace(0, 0, discreteSubdivisions)
    for i in range(0, int(discreteSubdivisions * lengthAToCenter / (lengthAToCenter + lengthCenterToC))) :
        y[i] = ra

    for i in range(int(discreteSubdivisions * lengthAToCenter / (lengthAToCenter + lengthCenterToC)), discreteSubdivisions):
        y[i] = -1 * rc

    y[499] = 0
    y[0] = 0
    plt.figure(figureNumber)
    plt.title(plotTitle)
    plt.ylabel("Shear Force [lbs]")
    plt.xlabel("x [in]")
    plt.plot(x, y)

    return [x, y]

def plot_moment_diagram(distance, shearForce, figureNumber, plotTitle):
    moment = np.linspace(0, 0, discreteSubdivisions)
    for i in range(0, discreteSubdivisions):
        moment[i] = np.trapz(shearForce[0 : i], distance[0 : i])
    plt.figure(figureNumber)
    plt.title(plotTitle)
    plt.ylabel("Moment [lbs-in]")
    plt.xlabel("x [in]")
    plt.plot(distance, moment)

    # Adjust by Kt
    for i in range(int(discreteSubdivisions * (lengthSec1 + lengthSec2) / (lengthAToCenter + lengthCenterToC)),
                   int(discreteSubdivisions * (lengthAToCenter + lengthSec3MiddleToEnd) / (lengthAToCenter + lengthCenterToC))):
        moment[i] = moment[i] * Kt
    plt.figure(figureNumber + 1)
    plt.title(plotTitle + " - Kt adjusted")
    plt.ylabel("Moment [lbs-in]")
    plt.xlabel("x [in]")
    plt.plot(distance, moment)
    return [distance, moment]

def cs_from_d(D): # From Table 5-4
    if D <= 0.3:
        return 1
    elif D <= 2 :
        return (D/0.3)**-0.11
    else:
        return 0.859 - 0.02125 * D


def get_section_diameter(torque, moment):
    return ((32 * N / math.pi) * ((Kt * moment / SnPrime) ** 2 + 0.75 * (torque / Sy)**2) ** 1/2) ** (1/3)

def get_section_diameters(tangentialMoment, radialMoment, torque):
    section1MaxMoment = 0
    for i in range(0, int(discreteSubdivisions * lengthSec1 / totalLength)):
        if math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2) > section1MaxMoment:
            section1MaxMoment = math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2)
    section1Diameter = get_section_diameter(torque, section1MaxMoment)

    section2MaxMoment = 0
    for i in range(int(discreteSubdivisions * lengthSec1 / totalLength),
                   int(discreteSubdivisions * (lengthSec1 + lengthSec2) / totalLength)):
        if math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2) > section2MaxMoment:
            section2MaxMoment = math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2)
    section2Diameter = get_section_diameter(torque, section2MaxMoment)

    section3MaxMoment = 0

    for i in range(int(discreteSubdivisions * (lengthSec1 + lengthSec2) / totalLength),
                   int(discreteSubdivisions * (lengthSec1 + lengthSec2 + lengthSec3MiddleToEnd + lengthSec3ToMiddle) / totalLength)):
        if math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2) > section3MaxMoment:
            section3MaxMoment = math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2)
    section3Diameter = get_section_diameter(torque, section3MaxMoment)

    section4MaxMoment = 0
    for i in range(int(discreteSubdivisions * (lengthSec1 + lengthSec2 + lengthSec3MiddleToEnd + lengthSec3ToMiddle) / totalLength),
                   int(discreteSubdivisions * (lengthSec1 + lengthSec2 + lengthSec3MiddleToEnd + lengthSec3ToMiddle + lengthSec4) / totalLength)):
        if math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2) > section4MaxMoment:
            section4MaxMoment = math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2)
    section4Diameter = get_section_diameter(0, section4MaxMoment)

    section5MaxMoment = 0
    for i in range(int(discreteSubdivisions * (lengthSec1 + lengthSec2 + lengthSec3MiddleToEnd + lengthSec3ToMiddle + lengthSec4) / totalLength),
                   int(discreteSubdivisions * (lengthSec1 + lengthSec2 + lengthSec3MiddleToEnd + lengthSec3ToMiddle + lengthSec4 + lengthSec5) / totalLength)):
        if math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2) > section5MaxMoment:
            section5MaxMoment = math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2)
    section5Diameter = get_section_diameter(0, section5MaxMoment)


    section6MaxMoment = 0
    for i in range(int(discreteSubdivisions * (lengthSec1 + lengthSec2 + lengthSec3MiddleToEnd + lengthSec3ToMiddle + lengthSec4 + lengthSec5) / totalLength),
                   int(discreteSubdivisions)):
        if math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2) > section6MaxMoment:
            section6MaxMoment = math.sqrt(tangentialMoment[i] ** 2 + radialMoment[i] ** 2)
    section6Diameter = get_section_diameter(0, section6MaxMoment)

    return [section1Diameter, section2Diameter, section3Diameter, section4Diameter, section5Diameter, section6Diameter]

inputSpeed = 1750
Kt = 3
Cr = 0.81 # 99% reliability from Table 5-3
D = 2
Sn = 310000 # SAE 9255 Q&T400 (Appendix 3) [ksi]
Sy = 287000 # SAE 9255 Q&T400 (Appendix 3) [ksi]
discreteSubdivisions = 500

Cs = cs_from_d(D)
SnPrime = Cs * Cr * Sn
N = 3 # Safety Factor

print("1) Section Lengths [in]: " + str([lengthSec1, lengthSec2, lengthSec3ToMiddle + lengthSec3MiddleToEnd, lengthSec4, lengthSec5, lengthSec6]))
print("2) Sn [psi]: " + str(Sn))
print("2) Sy [psi]: " + str(Sy))
print("3) Cs : " + str(Cs))
print("4) Sn' [psi]: " + str(SnPrime))
print("5) N: " + str(N))

angularVelocityIn = inputSpeed * 2 * math.pi / 60 # rpm
torque = powerIn / angularVelocityIn
print("6) Torque [lbs-in]: " + str(torque))

tangentialForceGear = torque / (pinionDiameter / 2)
radialForceGear = tangentialForceGear * math.tan(math.pi / 9)

print("7) Tangential Force Gear [lbs]: " + str(tangentialForceGear))
print("8) Radial Force Gear [lbs]: " + str(radialForceGear))

def father_of_get_section_diameters(tangentialForceGear, radialForceGear, startingFigureNumber, torque):
    [distance, shearForceTangential] = plot_shear_diagram(tangentialForceGear, lengthAToCenter, lengthCenterToC, startingFigureNumber, "Shaft 3 - Tangential Shear Force")
    [distance, shearForceRadial] = plot_shear_diagram(radialForceGear, lengthAToCenter, lengthCenterToC, startingFigureNumber + 1, "Shaft 3 - Radial Shear Force")

    [distance, tangentialMoment] = plot_moment_diagram(distance, shearForceTangential, startingFigureNumber + 2, "Shaft 3 - Moment Diagram Tangential")
    [distance, radialMoment] = plot_moment_diagram(distance, shearForceRadial, startingFigureNumber + 4, "Shaft 3 - Moment Diagram Radial")

    return get_section_diameters(tangentialMoment, radialMoment, torque)

section_diameters = father_of_get_section_diameters(tangentialForceGear, radialForceGear, 1, torque)
print("14) Section Diameters [in]: " + str(section_diameters))

print("\nSHAFT 3\n")

print("1) Section Lengths [in]: " + str([lengthSec1, lengthSec2, lengthSec3ToMiddle + lengthSec3MiddleToEnd, lengthSec4, lengthSec5, lengthSec6]))
print("2) Sn [psi]: " + str(Sn))
print("2) Sy [psi]: " + str(Sy))
print("3) Cs : " + str(Cs))
print("4) Sn' [psi]: " + str(SnPrime))
print("5) N: " + str(N))

angularVelocityIn = (inputSpeed / VR) * 2 * math.pi / 60
torque = powerIn / angularVelocityIn
print("6) Torque [lbs-in]: " + str(torque))

tangentialForceGear = torque / (gearDiameter / 2)
radialForceGear = tangentialForceGear * math.tan(math.pi / 9)
print("7) Tangential Force Gear [lbs]: " + str(tangentialForceGear))
print("8) Radial Force Gear [lbs]: " + str(radialForceGear))

section_diameters = father_of_get_section_diameters(tangentialForceGear, radialForceGear, 7, torque)
print("14) Output Section Diameters [in]: " + str(section_diameters))

print("\nSHAFT 2\n")

lengthA = 0.5
lengthB = 1
lengthC = 0.125
lengthD = 0.125
lengthEFirstHalf = 0.8565
lengthESecondHalf = 0.8565
lengthF = 6.074
lengthGFirstHalf = 0.8565
lengthGSecondHalf = 0.8565
lengthH = 0.125
lengthI = 0.125
lengthJ = 1
lengthK = 0.5
totalLength = lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf + lengthESecondHalf + lengthF + lengthGFirstHalf + \
    lengthGSecondHalf + lengthH + lengthI + lengthJ + lengthK

print("1) Section Lengths [in]: " + str([lengthA, lengthB, lengthC, lengthD, lengthEFirstHalf + lengthESecondHalf, lengthF,
                                         lengthGFirstHalf + lengthGSecondHalf, lengthH, lengthI, lengthJ, lengthK]))
print("2) Sn [psi]: " + str(Sn))
print("2) Sy [psi]: " + str(Sy))
print("3) Cs : " + str(Cs))
print("4) Sn' [psi]: " + str(SnPrime))
print("5) N: " + str(N))

shaftSpeed = (inputSpeed / VR) * 2 * math.pi / 60 # rpm
torque = power * 34000 / shaftSpeed
print("6) Torque [lbs-in]: " + str(torque))
tangentialForceInputGear = torque / (gearDiameter / 2)
tangentialForceOutputPinion = torque / (pinionDiameter / 2)

radialForceInputGear = tangentialForceGear * math.tan(math.pi / 9)
radialForceOutputPinion = tangentialForceOutputPinion * math.tan(math.pi / 9)

print("7) Tangential Force Input Gear [lbs]: " + str(tangentialForceInputGear))
print("7) Tangential Force Output Pinion [lbs]: " + str(tangentialForceOutputPinion))

print("8) Radial Force Input Gear [lbs]: " + str(radialForceInputGear))
print("8) Radial Force Output Pinion [lbs]: " + str(radialForceOutputPinion))

rkx = ((tangentialForceInputGear * (0.5 * lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf)) + \
      tangentialForceOutputPinion * (0.5 * lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf + lengthESecondHalf
                                     + lengthF + lengthGFirstHalf)) \
      / (totalLength - 0.5 * lengthA - 0.5 * lengthB)

rax = tangentialForceOutputPinion + tangentialForceInputGear - rkx

print("9) Tangential Reaction Force K: " + str(rkx))
print("9) Tangential Reaction Force A: " + str(rax))

rky = ((radialForceInputGear * (0.5 * lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf)) + \
      radialForceOutputPinion * (0.5 * lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf + lengthESecondHalf
                                     + lengthF + lengthGFirstHalf))\
      / (totalLength - 0.5 * lengthA - 0.5 * lengthB)
ray = radialForceOutputPinion + radialForceInputGear - rky

print("10) Radial Reaction Force K: " + str(rky))
print("10) Radial Reaction Force A: " + str(ray))

x = np.linspace(0, lengthAToCenter + lengthCenterToC, discreteSubdivisions)
shearForceTangential = np.linspace(0, 0, discreteSubdivisions)
for i in range(0, int(discreteSubdivisions * ((0.5 * lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf) / totalLength))):
    shearForceTangential[i] = rax

for i in range(int(discreteSubdivisions * ((0.5 * lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf) / totalLength)),
               int(discreteSubdivisions * (0.5 * lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf + lengthESecondHalf
                                     + lengthF + lengthGFirstHalf) / totalLength)):
    shearForceTangential[i] = rax - tangentialForceInputGear

for i in range(int(discreteSubdivisions * (0.5 * lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf + lengthESecondHalf
                                     + lengthF + lengthGFirstHalf) / totalLength), discreteSubdivisions):
    shearForceTangential[i] = rax - tangentialForceInputGear - tangentialForceOutputPinion

shearForceTangential[discreteSubdivisions - 1] = 0
shearForceTangential[0] = 0

plt.figure(13)
plt.title("Tangential Shear Force")
plt.ylabel("Shear Force [lbs]")
plt.xlabel("x [in]")
plt.plot(x, shearForceTangential)

shearForceRadial = np.linspace(0, 0, discreteSubdivisions)
for i in range(0, int(discreteSubdivisions * ((0.5 * lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf) / totalLength))):
    shearForceRadial[i] = ray

for i in range(int(discreteSubdivisions * ((0.5 * lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf) / totalLength)),
               int(discreteSubdivisions * (0.5 * lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf + lengthESecondHalf
                                     + lengthF + lengthGFirstHalf) / totalLength)):
    shearForceRadial[i] = ray - radialForceInputGear

for i in range(int(discreteSubdivisions * (0.5 * lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf + lengthESecondHalf
                                     + lengthF + lengthGFirstHalf) / totalLength), discreteSubdivisions):
    shearForceRadial[i] = ray - radialForceInputGear - radialForceOutputPinion

shearForceRadial[discreteSubdivisions - 1] = 0
shearForceRadial[0] = 0

plt.figure(14)
plt.title("Radial Shear Force")
plt.ylabel("Shear Force [lbs]")
plt.xlabel("x [in]")
plt.plot(x, shearForceRadial)

momentTangential = np.linspace(0, 0, discreteSubdivisions)
for i in range(0, discreteSubdivisions):
    momentTangential[i] = np.trapz(shearForceTangential[0: i], x[0: i])
plt.figure(15)
plt.title("Tangential Moment Diagram - Intermediate Shaft")
plt.ylabel("Moment [lbs-in]")
plt.xlabel("x [in]")
plt.plot(x, momentTangential)

# Adjust by Kt
for i in range(int(discreteSubdivisions * (lengthA + lengthB) / (totalLength)),
               int(discreteSubdivisions * (lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf + lengthESecondHalf) / (totalLength))):
    momentTangential[i] = momentTangential[i] * Kt

for i in range(int(discreteSubdivisions * (lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf + lengthESecondHalf + lengthF) / (totalLength)),
               int(discreteSubdivisions * (totalLength - lengthK - lengthJ) / (totalLength))):
    momentTangential[i] = momentTangential[i] * Kt

plt.figure(16)
plt.title("Tangential Moment Diagram - Intermediate Shaft - Kt adjusted")
plt.ylabel("Moment [lbs-in]")
plt.xlabel("x [in]")
plt.plot(x, momentTangential)

momentRadial = np.linspace(0, 0, discreteSubdivisions)
for i in range(0, discreteSubdivisions):
    momentRadial[i] = np.trapz(shearForceRadial[0: i], x[0: i])
plt.figure(17)
plt.title("Radial Moment Diagram - Intermediate Shaft")
plt.ylabel("Moment [lbs-in]")
plt.xlabel("x [in]")
plt.plot(x, momentRadial)

# Adjust by Kt
for i in range(int(discreteSubdivisions * (lengthA + lengthB) / (totalLength)),
               int(discreteSubdivisions * (lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf + lengthESecondHalf) / (totalLength))):
    momentRadial[i] = momentRadial[i] * Kt

for i in range(int(discreteSubdivisions * (lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf + lengthESecondHalf + lengthF) / (totalLength)),
               int(discreteSubdivisions * (totalLength - lengthK - lengthJ) / (totalLength))):
    momentRadial[i] = momentRadial[i] * Kt

plt.figure(18)
plt.title("Radial Moment Diagram - Intermediate Shaft - Kt adjusted")
plt.ylabel("Moment [lbs-in]")
plt.xlabel("x [in]")
plt.plot(x, momentRadial)

sectionAMaxMoment = 0
for i in range(0, int(discreteSubdivisions * lengthA / totalLength)):
    if math.sqrt(momentTangential[i] ** 2 + momentRadial[i] ** 2) > sectionAMaxMoment:
        sectionAMaxMoment = math.sqrt(momentTangential[i] ** 2 + momentRadial[i] ** 2)
sectionADiameter = get_section_diameter(0, sectionAMaxMoment)

sectionBMaxMoment = 0
for i in range(int(discreteSubdivisions * lengthA / totalLength),
               int(discreteSubdivisions * (lengthA + lengthB) / totalLength)):
    if math.sqrt(momentTangential[i] ** 2 + momentRadial[i] ** 2) > sectionBMaxMoment:
        sectionBMaxMoment = math.sqrt(momentTangential[i] ** 2 + momentRadial[i] ** 2)
sectionBDiameter = get_section_diameter(0, sectionBMaxMoment)

sectionCDEMaxMoment = 0
for i in range(int(discreteSubdivisions * (lengthA + lengthB) / totalLength),
               int(discreteSubdivisions * (lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf + lengthESecondHalf) / totalLength)):
    if math.sqrt(momentTangential[i] ** 2 + momentRadial[i] ** 2) > sectionCDEMaxMoment:
        sectionCDEMaxMoment = math.sqrt(momentTangential[i] ** 2 + momentRadial[i] ** 2)
sectionCDEDiameter = get_section_diameter(torque, sectionCDEMaxMoment)

sectionFMaxMoment = 0
for i in range(int(discreteSubdivisions * (lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf + lengthESecondHalf)),
               int(discreteSubdivisions * (lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf + lengthESecondHalf + lengthF) / totalLength)):
    if math.sqrt(momentTangential[i] ** 2 + momentRadial[i] ** 2) > sectionFMaxMoment:
        sectionFMaxMoment = math.sqrt(momentTangential[i] ** 2 + momentRadial[i] ** 2)
sectionFDiameter = get_section_diameter(torque, sectionFMaxMoment)

sectionGHIMaxMoment = 0
for i in range(int(discreteSubdivisions * (lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf + lengthESecondHalf + lengthF) / totalLength),
               int(discreteSubdivisions * (lengthA + lengthB + lengthC + lengthD + lengthEFirstHalf + lengthESecondHalf + lengthF + lengthGFirstHalf + lengthGSecondHalf + lengthH + lengthI) / totalLength)):
    if math.sqrt(momentTangential[i] ** 2 + momentRadial[i] ** 2) > sectionGHIMaxMoment:
        sectionGHIMaxMoment = math.sqrt(momentTangential[i] ** 2 + momentRadial[i] ** 2)
sectionGHIDiameter = get_section_diameter(torque, sectionGHIMaxMoment)

sectionJMaxMoment = 0
for i in range(int(discreteSubdivisions * (totalLength - lengthK - lengthJ) / totalLength),
               int(discreteSubdivisions * (totalLength - lengthK) / totalLength)):
    if math.sqrt(momentTangential[i] ** 2 + momentRadial[i] ** 2) > sectionJMaxMoment:
        sectionJMaxMoment = math.sqrt(momentTangential[i] ** 2 + momentRadial[i] ** 2)
sectionJDiameter = get_section_diameter(0, sectionJMaxMoment)

sectionKMaxMoment = 0
for i in range(int(discreteSubdivisions * (totalLength - lengthK) / totalLength),
               int(discreteSubdivisions * (totalLength) / totalLength)):
    if math.sqrt(momentTangential[i] ** 2 + momentRadial[i] ** 2) > sectionKMaxMoment:
        sectionKMaxMoment = math.sqrt(momentTangential[i] ** 2 + momentRadial[i] ** 2)
sectionKDiameter = get_section_diameter(0, sectionKMaxMoment)

intermediateShaftOutputSection = [sectionADiameter, sectionBDiameter, sectionCDEDiameter, sectionFDiameter, sectionGHIDiameter, sectionJDiameter, sectionKDiameter]
print("14) Intermediate Shaft Diameters [in]: " + str(intermediateShaftOutputSection))

plt.show()
