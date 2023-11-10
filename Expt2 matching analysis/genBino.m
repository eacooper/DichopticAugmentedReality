%given the weight assigned to the higher contrast eye and the left/right eye's contrast, find the predicted
%binocular contrast
function predicts = genBino(lowC,highC,weights)
%weights is a list so that each column of "predicts" will contain the predicted value to one weight value
predicts = [];

for w = weights
    HighW = w; %high contrast weight
    LowW = 1-w; %low contrast weight
    Model_resp = HighW.*highC + LowW.*lowC; %weighted combination
    predicts = [predicts, Model_resp]; %add this weight's prediction as a new column
end

