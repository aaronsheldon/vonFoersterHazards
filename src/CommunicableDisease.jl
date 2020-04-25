"""
    CommunicableDisease

Demonstration implementation of a minimum module of age dependent 
communicable diseases outlined in the notebook. This module effectively
infers the counterfactual utilization in the asymptotic limit of infinite
hospital capacity.
"""
module CommunicableDisease

using vonFoersterHazards

function gompertzrate end
function makehamrate end
function recoveryrate end
function hospitalizationrate end
function ageingmortality end
function infectionmortality end
function hazardrate end
function hazardrate end
function scatterrate end
function birthrate end
function runengine end

end