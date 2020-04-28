"""
    CommunicableDisease

Demonstration implementation of a minimum module of age dependent 
communicable diseases outlined in the notebook. This module can be
used to infer the counterfactual utilization under a variety of
hospital capacities. Remember that our transition matrix convention
is that SOURCES are columns and DESTINATIONS are rows.
"""
module CommunicableDisease

using vonFoersterHazards

export stateontology,
       statetransitions,
       parameterontology,
       parameterestimates

const stateontology = [ "Person"    "Susceptible";
                        "Person"    "Infected";
                        "Person"    "Hospitalized";
                        "Person"    "Recovered";
                        "Person"    "Discharged";
                        "Person"    "Ageing Fatalities";
                        "Person"    "Infected Fatalities";
                        "Person"    "Hospitalized Fatalities";
                        "Acute Bed" "Available";
                        "Acute Bed" "Occupied";
                        "Acute Bed" "Refractory" ]

const statetransitions = [  1  2;
                            1  6;
                            2  3;
                            2  4;
                            2  7;
                            3  5;
                            3  8;
                            4  6;
                            5  6;
                            9 10;
                           10 11;
                           11  9 ]

const parameterontology = [ "Endogenous" "Background Mortality" "Base Rate"        "Deaths per Person Year";
                            "Endogenous" "Background Mortality" "Doubling Time"    "per Age Year";
                            "Endogenous" "Background Admission" "Base Rate"        "Admissions per Person Year";
                            "Endogenous" "Background Admission" "Doubling Time"    "per Age Year";
                            "Endogenous" "Background Discharge" "Base Rate"        "Discharges per Person Day";
                            "Endogenous" "Background Discharge" "Half Life"        "per Age Year";
                            "Endogenous" "Puberty"              "Doubling Time"    "per Age Year";
                            "Endogenous" "Puberty"              "Mid-Point"        "Age Year";
                            "Endogenous" "Infection"            "Acceleration"     "Dimensionless";
                            "Exogenous"  "Infection"            "Contact Rate"     "Contact per Person Day";
                            "Exogenous"  "Infection"            "Supression Ratio" "Dimensionless";
                            "Endogenous" "Acute Bed"            "Refraction Rate"  "Availabilities per Bed Day"; ]

const parameterestimates = [ 1.0/32000.0;
                                     7.0;
                               1.0/320.0;
                                    14.0;
                                1.0/16.0;
                                    49.0; 
                                     2.0;
                                    16.0;
                                     1.3;
                                     0.1;
                                     0.3;
                                     2.0; ]

function gompertzrate end
function makehamrate end
function recoveryrate end
function hospitalizationrate end
function ageingmortality end
function infectionmortality end
function hazardrate end
function scatterrate end
function birthrate end
function runengine end

end