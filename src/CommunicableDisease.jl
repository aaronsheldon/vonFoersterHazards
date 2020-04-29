"""
    CommunicableDisease

Demonstration implementation of a minimum module of age dependent 
communicable diseases outlined in the notebook. This module can be
used to infer the counterfactual utilization under a variety of
hospital capacities. Remember that our transition matrix convention
is that SOURCES are columns and TRAGETS are rows.
"""
module CommunicableDisease

using vonFoersterHazards

export stateontology,
       statetransitions,
       parameterontology,
       parameterestimates,
       hazardrate,
       scatterrate,
       birthrate,
       runengine

const stateontology = ( head = [ "Noun"      "Adjective" ],
                        body = [ "Person"    "Susceptible";
                                 "Person"    "Infected";
                                 "Person"    "Hospitalized";
                                 "Person"    "Recovered";
                                 "Person"    "Discharged";
                                 "Person"    "Ageing Fatality";
                                 "Person"    "Infected Fatality";
                                 "Person"    "Hospitalized Fatality";
                                 "Acute Bed" "Available";
                                 "Acute Bed" "Occupied";
                                 "Acute Bed" "Refractory" ] )

const statetransitions = ( head = [ "Source" "Target" ],
                           body = [        1        2;
                                           1        6;
                                           2        3;
                                           2        4;
                                           2        7;
                                           3        5;
                                           3        8;
                                           4        6;
                                           5        6;
                                           9       10;
                                          10       11;
                                          11        9 ] )

const parameterontology = ( head = [ "Decomposition" "Process"              "Metric"           "Units" ],
                            body = [ "Endogenous"    "Background Mortality" "Base Rate"        "Deaths per Person Year";
                                     "Endogenous"    "Background Mortality" "Doubling Time"    "per Age Year";
                                     "Endogenous"    "Background Admission" "Base Rate"        "Admissions per Person Year";
                                     "Endogenous"    "Background Admission" "Doubling Time"    "per Age Year";
                                     "Endogenous"    "Background Discharge" "Base Rate"        "Discharges per Person Day";
                                     "Endogenous"    "Background Discharge" "Half Life"        "per Age Year";
                                     "Endogenous"    "Puberty"              "Doubling Time"    "per Age Year";
                                     "Endogenous"    "Puberty"              "Mid-Point"        "Age Year";
                                     "Endogenous"    "Infection"            "Acceleration"     "Dimensionless";
                                     "Exogenous"     "Infection"            "Contact Rate"     "Contacts per Person Day";
                                     "Exogenous"     "Infection"            "Supression Ratio" "Dimensionless";
                                     "Endogenous"    "Acute Bed"            "Refraction Rate"  "Availabilities per Bed Day"; ] )

const parameterestimates = ( head = [ "Value" ], 
                             body = [ 1.0/32000.0;
                                              7.0;
                                        1.0/320.0;
                                             14.0;
                                         1.0/16.0;
                                             49.0; 
                                              2.0;
                                             16.0;
                                              1.3;
                                         1.0/10.0;
                                          1.0/3.0;
                                              2.0; ] )
"""
    markdowntable(t)

Throw away function to generate a markdown table from the constants
that describe the model.
"""
function markdowntable(t)
    s = join(["|" * string(h) for h in t.head]) * "|"
    s *= "\n" * join(["|-" for h in t.head]) * "|"
    for r in eachindex(t.body[:, 1])
        s *= "\n" * join(["|" * string(b) for b in t.body[r, :]]) * "|"
    end
    display("text/markdown", s)
end

function hazardrate end
function scatterrate end
function birthrate end
function runengine end

end