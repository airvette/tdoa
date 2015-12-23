function locSource = MULocate(evVal)
% Takes in a 3x3 array that contains differing arrival times and the
%  locations of three WADN nodes in order to approximate the location of an
%  unknown radiation source (i.e. a malicious user).  
% The Time Difference of Arrival (TDOA) Data Fusion method from Sayed, 2005
%  is used to find the location of the user.  Some of the terminology used 
%  is as follows:
%  - x1, x2, x3: the x coordinate of nodes 1, 2 or 3 (x1 is set to 0)
%  - y1, y2, y3: the y coordinate of nodes 1, 2 or 3 (y1 is set to 0)
%  - delT21, delT31: the time difference in signal arrival between nodes 2
%    and 1 and nodes 3 and 1
%  - c: the speed of light, in m/us
%  - r21, r31: the range differences from the unknown source between nodes
%    2 and 1 and 3 and 1
%================ Inputs =========================================
% evVal = A 3x3 array that has the values for the most recent single
%  event that three or more WADN nodes reported on.  These values are
%  needed to determine the non-WADN source (Malicious User) location.  Each
%  row is identical and each row has the format:
%  [Unknown Sig St Time, Node X-Location, Node Y-Location]
%================ Outputs ========================================
% xSource = x value of the radiation source
% ySource = y value of the radiation source
%////////////////////////////////////////////////////////////////
% Jeffrey Guido, UCCS Masters Thesis, FA 2013
 
if min(evVal(:,1)) > 0 % Check first if valid entries are in the input
% if max(evVal(:,1)) > 0 % Check first if valid entries are in the input
    % Determine which node is closest to the source, S
    [~,first] = min(evVal(:,1));
    % Determine the order of the other nodes
    switch first 
        case 1
            [~, second] = min([inf evVal(2,1) evVal(3,1)]);
            [~, third] = max([-inf evVal(2,1) evVal(3,1)]);
            if second == third % if second equals third, then it is both the 
                % max and min and evVal(2,1) is equal to evVal(3,1)
                third = 3; % so set it to the element that max() gets to second
            end
        case 2
            [~, second] = min([evVal(1,1) inf evVal(3,1)]);
            [~, third] = max([evVal(1,1) -inf evVal(3,1)]);
            if second == third % if second equals third, then it is both the 
                % max and min and evVal(2,1) is equal to evVal(3,1)
                third = 3; % so set it to the element that max() gets to second
            end
        otherwise
            [~, second] = min([evVal(1,1) evVal(2,1) inf]);
            [~, third] = max([evVal(1,1) evVal(2,1) -inf]);
            if second == third % if second equals third, then it is both the 
                % max and min and evVal(2,1) is equal to evVal(3,1)
                third = 2; % so set it to the element that max() gets to second
            end
    end
 
    % To reduce the number of unknown varibles in the least squares equation
    %  below, the origin of the coordinate system is shifted to the location of 
    %  the first node to receive the signal from the source.  All other
    %  node coordinates are updated to match this shift
    x2 = evVal(second,2) - evVal(first,2); % Shift x2 location
    y2 = evVal(second,3) - evVal(first,3); % Shift y2 location
    x3 = evVal(third,2) - evVal(first,2); % Shift x3 location
    y3 = evVal(third,3) - evVal(first,3); % Shift y3 location
    % Determine time of arrival differences
    delT21 = evVal(second,1) - evVal(first,1); % Difference between 2 and 1
    delT31 = evVal(third,1) - evVal(first,1); % Difference between 3 and 1
 
    c = 300; % Speed of light in m/us
 
    % Determine the differences in range
    r21 = delT21 * c; % Range difference between nodes 2 and 1
    r31 = delT31 * c; % Range difference between nodes 3 and 1
 
    H = [x2 y2; x3 y3]; % Matrix H from Say05, pg 4, eqn 7
    C = [-r21; -r31]; % Vector C from Say05, pg 5, eqn 11
    K2sq = x2^2 + y2^2; % K2 Squared from Say05, pg4, eqn 6
    K3sq = x3^2 + y3^2; % K3 Squared from Say05, pg4, eqn 6
    D = .5.*[K2sq - r21^2; K3sq - r31^2]; % Vector D from Say05, pg 5, eqn 11
 
    Xint = [H\C, H\D]; % Intermediate values for vector X.  X 
        % contains the values [xm = xa*r1 + xb; ym = yar1 + yb] in the form of
        % [xa xb; ya yb] which is the location of the source when r1 is solved
 
    % The above expression, when substituted back into Say05, eqn 2, which is 
    %  r1^2 = xm^2 + ym^2, yields a polynomial of power 2.  The following line
    %  reorganizes Xint above into that polynomial so that the roots can be
    %  determined
    p = [(Xint(1,1)^2 + Xint(2,1)^2 - 1), ... % coefficient for r1^2
         (Xint(1,1)*Xint(1,2)*2 + Xint(2,1)*Xint(2,2)*2), ... % coeff for r1^1
         (Xint(1,2)^2 + Xint(2,2)^2)]; % coefficient for r1^0
    root = roots(p); % determine the roots of the polynomial p
 
    % Pick the largest root to be r1
    r1 = max(real(root));
    xm = r1*Xint(1,1) + Xint(1,2); % Calculate the location of xm in the 
        % shifted coordinate system
    ym = r1*Xint(2,1) + Xint(2,2);% Calculate the location of xm in the 
        % shifted coordinate system
    xSource = xm + evVal(first,2); % Shift x coord back to actual location
    ySource = ym + evVal(first,3); % Shift x coord back to actual location
 
    locSource = [xSource; ySource]; % Enter the source locations into a vector
else
    locSource = [-1;-1]; % Default output.  An invalid location since no 
        % nodes exist in the negative regions
end 
Status API Training Shop Blog About Pricing
© 2015 GitHub, Inc. Terms Privacy Security Contact Help
