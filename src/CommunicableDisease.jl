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

const ontology = [ "Person"    "Susceptible";
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

const transitions = [  1  2;
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