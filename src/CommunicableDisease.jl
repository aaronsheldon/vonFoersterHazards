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
       processontology,
       statetransitions,
       parameterontology,
       parameterestimates,
       configurationontology,
       configurationvalues,
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
                                 "Person"    "Hospital Fatality";
                                 "Acute Bed" "Available";
                                 "Acute Bed" "Occupied";
                                 "Acute Bed" "Refractory" ] )

const processontology = ( head = [ "Process"              "Source Noun" "Source Adjective" "Target Noun" "Target Adjective"],
                          body = [ "Background Mortality" "Person"      "Susceptible"      "Person"      "Ageing Fatality"; 
                                   "Background Mortality" "Person"      "Recovered"        "Person"      "Ageing Fatality";
                                   "Background Mortality" "Person"      "Discharged"       "Person"      "Ageing Fatality";
                                   "Infected Mortality"   "Person"      "Infected"         "Person"      "Infected Fatality"; 
                                   "Infected Mortality"   "Person"      "Hospitalized"     "Person"      "Hospital Fatality"; 
                                   "Admission"            "Person"      "Infected"         "Person"      "Hospitalized"; 
                                   "Discharge & Recovery" "Person"      "Infected"         "Person"      "Recovered"; 
                                   "Discharge & Recovery" "Person"      "Hospitalized"     "Person"      "Discharged";
                                   "Infection"            "Person"      "Susceptible"      "Person"      "Infected"; 
                                   "Bed Allocate"         "Acute Bed"   "Available"        "Acute Bed"   "Occupied"; 
                                   "Bed Release"          "Acute Bed"   "Occupied"         "Acute Bed"   "Refractory";
                                   "Bed Recover"          "Acute Bed"   "Refractory"       "Acute Bed"   "Available" ] )

const statetransitions = ( head = [ "Source" "Target" "Process" ],
                           body = [        1        2         5;
                                           1        6         1;
                                           2        3         3;
                                           2        4         4;
                                           2        7         2;
                                           3        5         4;
                                           3        8         2;
                                           4        6         1;
                                           5        6         1;
                                           9       10         6;
                                          10       11         7;
                                          11        9         8 ] )

const parameterontology = ( head = [ "Decomposition" "Process"              "Metric"            "Units" ],
                            body = [ "Endogenous"    "Background Mortality" "Base Rate"         "Deaths per Person Year";
                                     "Endogenous"    "Background Mortality" "Doubling Time"     "per Age Year";
                                     "Endogenous"    "Background Admission" "Base Rate"         "Admissions per Person Year";
                                     "Endogenous"    "Background Admission" "Doubling Time"     "per Age Year";
                                     "Endogenous"    "Background Discharge" "Base Rate"         "Discharges per Person Day";
                                     "Endogenous"    "Background Discharge" "Half Life"         "per Age Year";
                                     "Endogenous"    "Puberty"              "Doubling Time"     "per Age Year";
                                     "Endogenous"    "Puberty"              "Mid-Point"         "Age Year";
                                     "Endogenous"    "Infection"            "Acceleration"      "Dimensionless";
                                     "Exogenous"     "Infection"            "Contact Rate"      "Contacts per Person Day";
                                     "Exogenous"     "Infection"            "Suppression Ratio" "Dimensionless";
                                     "Endogenous"    "Acute Bed"            "Refraction Rate"   "Availabilities per Bed Day" ] )

const parameterestimates = ( head = [     "Value" ], 
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
                                              2.0 ] )

const configurationontology = ( head = [ "Name"             "Description"                        "Units" ], 
                                body = [ "Time Step"        "Eapsed time of a single iteration." "Days";
                                         "Step Count"       "Total iterations to evolve."        "Dimensionless";
                                         "Cohort Gestation" "Elapsed time to add a new cohort."  "Days" ] )

const configurationvalues = ( head = [ "Value" ],
                              body = [       1;
                                           730;
                                           365 ] )

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

"""
    hazardrate(p)

The hazard rate matrix is spectrally decomposed into 8 age-Eigen function matrices,
each mapping 11 states to 11 states.
"""
function hazardrate(p::population{Int64, Vector{Int64}, Matrix{Int64}, Array{Float64, 3}})
    zeros(Float64, 8, 11,11)
end

"""
    scatterrate(a)

The scattering rate decomposes the hazard rate into 8 age scattering cross sections.
"""
function scatterrate(a::Int64)
    zeros(Float64, 8)
end

"""
    birthrate(p)

The birth rate is constant and zero for all states except the first, susceptible, which
contains the daily average births.
"""
function birthrate(p::population{Int64, Vector{Int64}, Matrix{Int64}, Array{Float64, 3}})
    zeros(Int64, 11)
end

"""
    runengine

Run the simulation. At each step coerce the returned population into a dimensional
vector store and write to disk.
"""
function runengine end

end