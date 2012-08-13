function [x,w] = kpn_knots(n)

% [x,w] = kpn_knots(n)
%
% returns knots and weights for kpn for given number of points


% first recover the i-level
i=knots2lev_kpn(n);


if isnan(i)
   error(strcat('this number of points is not available:',num2str(n)))
else
    % now access the knots and weights table using the i2l map
    [x_t,w_t]=kpn_tabulated(kpn_lev2l_map(i));
    x = redistribute_knots(x_t);
    w = redistribute_weights(w_t);
end

% sort knots increasingly and weights accordingly
[x,sorter]=sort(x);
w=w(sorter);

end

%--------------------------------------------------------------------------------------------------------------------------------
% secondary functions
%--------------------------------------------------------------------------------------------------------------------------------


%----------------------------------------------------------------
function i = knots2lev_kpn(nb_knots)

% i = knots2lev_kpn(nb_knots)
%
% given the number of points needed, returns the corresponding condensed level i,
% according to kpn_lev_table.

switch nb_knots
    case 0
        i=kpn_lev_table(1,1);
    case 1
        i=kpn_lev_table(2,1);
    case 3
        i=kpn_lev_table(3,1);
    case 9
        i=kpn_lev_table(4,1);
    case 19
        i=kpn_lev_table(5,1);
    case 35
        i=kpn_lev_table(6,1);
    otherwise
        i=NaN;
end

end



%----------------------------------------------------------------
function l = kpn_lev2l_map(i)

% l = kpn_lev2l_map(i)
%
% returns the l-level (``true level'') corresponding to the level i (``condensed level''),
% as from table in kpn_lev_table

if i>5
    error('this level is not tabulated')
else
    % i to l map
    l=kpn_lev_table(i+1,2);
    return 
end
end


%----------------------------------------------------------------
function out = redistribute_knots(in)

nin=length(in);
nout=1+2*(nin-1);

out=zeros(1,nout);

out(1)=in(1);

in(1)=[];

T=[-in;in];

out(2:end)=T(:);

end


%----------------------------------------------------------------
function out = redistribute_weights(in)

nin=length(in);
nout=1+2*(nin-1);

out=zeros(1,nout);

out(1)=in(1);

in(1)=[];

T=[in;in];

out(2:end)=T(:);

end



%----------------------------------------------------------------
function [ n, w ] = kpn_tabulated ( l )

% [n,w] = kpn_tabulated ( l )
%
% returns knots and weights kpn at level l from tabulated values.
% l has to be in the interval 1 <= l <= 25
% 
%
% type edit kpn_tabulated for credits and references



% KPN provides data for Kronrod-Patterson quadrature with a normal weight.
%
%  Discussion:
%
%    This data assumes integration over the interval (-oo,+oo) with 
%    weight function w(x) = exp(-x*x/2)/sqrt(2*pi).
%
%    For all orders L, the rule is formed by
%      X(1) with weight W(1),
%      X(2) and -X(2) with weight W(2),
%      X(3) and -X(3) with weight W(3) and so on.
%
%  Modified:
%
%    01 June 2010
%
%  Author:
%
%    Florian Heiss, Viktor Winschel
%
%  Reference:
%
%    Florian Heiss, Viktor Winschel,
%    Likelihood approximation by numerical integration on sparse grids,
%    Journal of Econometrics,
%    Volume 144, 2008, pages 62-80.
%
%    Alan Genz, Bradley Keister,
%    Fully symmetric interpolatory rules for multiple integrals
%    over infinite regions with Gaussian weight,
%    Journal of Computational and Applied Mathematics,
%    Volume 71, 1996, pages 299-309.
%
%    Thomas Patterson,
%    The optimal addition of points to quadrature formulae,
%    Mathematics of Computation,
%    Volume 22, Number 104, October 1968, pages 847-856.
%
%  Parameters:
%
%    Input, integer L, the level of the rule.
%
%    Output, real N(L), the nodes.
%
%    Output, real W(L), the weights.
%
  switch l
  case 1
    n = [0.0000000000000000e+000];
    w = [1.0000000000000000e+000];
  case 2
    n = [0.0000000000000000e+000; 1.7320508075688772e+000];
    w = [6.6666666666666663e-001; 1.6666666666666666e-001];
  case 3
    n = [0.0000000000000000e+000; 1.7320508075688772e+000];
    w = [6.6666666666666674e-001; 1.6666666666666666e-001];
  case 4
    n = [0.0000000000000000e+000; 7.4109534999454085e-001; 1.7320508075688772e+000; 4.1849560176727323e+000];
    w = [4.5874486825749189e-001; 1.3137860698313561e-001; 1.3855327472974924e-001; 6.9568415836913987e-004];
  case 5
    n = [0.0000000000000000e+000; 7.4109534999454085e-001; 1.7320508075688772e+000; 2.8612795760570582e+000; 4.1849560176727323e+000];
    w = [2.5396825396825407e-001; 2.7007432957793776e-001; 9.4850948509485125e-002; 7.9963254708935293e-003; 9.4269457556517470e-005];
  case 6
    n = [0.0000000000000000e+000; 7.4109534999454085e-001; 1.7320508075688772e+000; 2.8612795760570582e+000; 4.1849560176727323e+000];
    w = [2.5396825396825429e-001; 2.7007432957793776e-001; 9.4850948509485070e-002; 7.9963254708935293e-003; 9.4269457556517551e-005];
  case 7
    n = [0.0000000000000000e+000; 7.4109534999454085e-001; 1.7320508075688772e+000; 2.8612795760570582e+000; 4.1849560176727323e+000];
    w = [2.5396825396825418e-001; 2.7007432957793781e-001; 9.4850948509485014e-002; 7.9963254708935311e-003; 9.4269457556517592e-005];
  case 8
    n = [0.0000000000000000e+000; 7.4109534999454085e-001; 1.7320508075688772e+000; 2.8612795760570582e+000; 4.1849560176727323e+000];
    w = [2.5396825396825418e-001; 2.7007432957793781e-001; 9.4850948509485042e-002; 7.9963254708935276e-003; 9.4269457556517375e-005];
  case 9
    n = [0.0000000000000000e+000; 7.4109534999454085e-001; 1.2304236340273060e+000; 1.7320508075688772e+000; 2.5960831150492023e+000; 2.8612795760570582e+000; 4.1849560176727323e+000; 5.1870160399136562e+000; 6.3633944943363696e+000];
    w = [2.6692223033505302e-001; 2.5456123204171222e-001; 1.4192654826449365e-002; 8.8681002152028010e-002; 1.9656770938777492e-003; 7.0334802378279075e-003; 1.0563783615416941e-004; -8.2049207541509217e-007; 2.1136499505424257e-008];
  case 10
    n = [0.0000000000000000e+000; 7.4109534999454085e-001; 1.2304236340273060e+000; 1.7320508075688772e+000; 2.5960831150492023e+000; 2.8612795760570582e+000; 3.2053337944991944e+000; 4.1849560176727323e+000; 5.1870160399136562e+000; 6.3633944943363696e+000];
    w = [3.0346719985420623e-001; 2.0832499164960877e-001; 6.1151730125247716e-002; 6.4096054686807610e-002; 1.8085234254798462e-002; -6.3372247933737571e-003; 2.8848804365067559e-003; 6.0123369459847997e-005; 6.0948087314689840e-007; 8.6296846022298632e-010];
  case 11
    n = [0.0000000000000000e+000; 7.4109534999454085e-001; 1.2304236340273060e+000; 1.7320508075688772e+000; 2.5960831150492023e+000; 2.8612795760570582e+000; 3.2053337944991944e+000; 4.1849560176727323e+000; 5.1870160399136562e+000; 6.3633944943363696e+000];
    w = [3.0346719985420623e-001; 2.0832499164960872e-001; 6.1151730125247709e-002; 6.4096054686807541e-002; 1.8085234254798459e-002; -6.3372247933737545e-003; 2.8848804365067555e-003; 6.0123369459847922e-005; 6.0948087314689830e-007; 8.6296846022298839e-010];
  case 12
    n = [0.0000000000000000e+000; 7.4109534999454085e-001; 1.2304236340273060e+000; 1.7320508075688772e+000; 2.5960831150492023e+000; 2.8612795760570582e+000; 3.2053337944991944e+000; 4.1849560176727323e+000; 5.1870160399136562e+000; 6.3633944943363696e+000];
    w = [3.0346719985420623e-001; 2.0832499164960872e-001; 6.1151730125247716e-002; 6.4096054686807624e-002; 1.8085234254798466e-002; -6.3372247933737545e-003; 2.8848804365067559e-003; 6.0123369459847841e-005; 6.0948087314689830e-007; 8.6296846022298963e-010];
  case 13
    n = [0.0000000000000000e+000; 7.4109534999454085e-001; 1.2304236340273060e+000; 1.7320508075688772e+000; 2.5960831150492023e+000; 2.8612795760570582e+000; 3.2053337944991944e+000; 4.1849560176727323e+000; 5.1870160399136562e+000; 6.3633944943363696e+000];
    w = [3.0346719985420600e-001; 2.0832499164960883e-001; 6.1151730125247730e-002; 6.4096054686807638e-002; 1.8085234254798459e-002; -6.3372247933737580e-003; 2.8848804365067555e-003; 6.0123369459847868e-005; 6.0948087314689830e-007; 8.6296846022298756e-010];
  case 14
    n = [0.0000000000000000e+000; 7.4109534999454085e-001; 1.2304236340273060e+000; 1.7320508075688772e+000; 2.5960831150492023e+000; 2.8612795760570582e+000; 3.2053337944991944e+000; 4.1849560176727323e+000; 5.1870160399136562e+000; 6.3633944943363696e+000];
    w = [3.0346719985420617e-001; 2.0832499164960874e-001; 6.1151730125247702e-002; 6.4096054686807596e-002; 1.8085234254798459e-002; -6.3372247933737563e-003; 2.8848804365067555e-003; 6.0123369459847936e-005; 6.0948087314689851e-007; 8.6296846022298322e-010];
  case 15
    n = [0.0000000000000000e+000; 7.4109534999454085e-001; 1.2304236340273060e+000; 1.7320508075688772e+000; 2.5960831150492023e+000; 2.8612795760570582e+000; 3.2053337944991944e+000; 4.1849560176727323e+000; 5.1870160399136562e+000; 6.3633944943363696e+000];
    w = [3.0346719985420612e-001; 2.0832499164960874e-001; 6.1151730125247723e-002; 6.4096054686807652e-002; 1.8085234254798459e-002; -6.3372247933737597e-003; 2.8848804365067563e-003; 6.0123369459848091e-005; 6.0948087314689851e-007; 8.6296846022298983e-010];
  case 16
    n = [0.0000000000000000e+000; 2.4899229757996061e-001; 7.4109534999454085e-001; 1.2304236340273060e+000; 1.7320508075688772e+000; 2.2336260616769419e+000; 2.5960831150492023e+000; 2.8612795760570582e+000; 3.2053337944991944e+000; 3.6353185190372783e+000; 4.1849560176727323e+000; 5.1870160399136562e+000; 6.3633944943363696e+000; 7.1221067008046166e+000; 7.9807717985905606e+000; 9.0169397898903032e+000];
    w = [2.5890005324151566e-001; 2.8128101540033167e-002; 1.9968863511734550e-001; 6.5417392836092561e-002; 6.1718532565867179e-002; 1.7608475581318002e-003; 1.6592492698936010e-002; -5.5610063068358157e-003; 2.7298430467334002e-003; 1.5044205390914219e-005; 5.9474961163931621e-005; 6.1435843232617913e-007; 7.9298267864869338e-010; 5.1158053105504208e-012; -1.4840835740298868e-013; 1.2618464280815118e-015];
  case 17
    n = [0.0000000000000000e+000; 2.4899229757996061e-001; 7.4109534999454085e-001; 1.2304236340273060e+000; 1.7320508075688772e+000; 2.2336260616769419e+000; 2.5960831150492023e+000; 2.8612795760570582e+000; 3.2053337944991944e+000; 3.6353185190372783e+000; 4.1849560176727323e+000; 5.1870160399136562e+000; 5.6981777684881099e+000; 6.3633944943363696e+000; 7.1221067008046166e+000; 7.9807717985905606e+000; 9.0169397898903032e+000];
    w = [1.3911022236338039e-001; 1.0387687125574284e-001; 1.7607598741571459e-001; 7.7443602746299481e-002; 5.4677556143463042e-002; 7.3530110204955076e-003; 1.1529247065398790e-002; -2.7712189007789243e-003; 2.1202259559596325e-003; 8.3236045295766745e-005; 5.5691158981081479e-005; 6.9086261179113738e-007; -1.3486017348542930e-008; 1.5542195992782658e-009; -1.9341305000880955e-011; 2.6640625166231651e-013; -9.9313913286822465e-016];
  case 18
    n = [0.0000000000000000e+000; 2.4899229757996061e-001; 7.4109534999454085e-001; 1.2304236340273060e+000; 1.7320508075688772e+000; 2.2336260616769419e+000; 2.5960831150492023e+000; 2.8612795760570582e+000; 3.2053337944991944e+000; 3.6353185190372783e+000; 4.1849560176727323e+000; 4.7364330859522967e+000; 5.1870160399136562e+000; 5.6981777684881099e+000; 6.3633944943363696e+000; 7.1221067008046166e+000; 7.9807717985905606e+000; 9.0169397898903032e+000];
    w = [5.1489450806921377e-004; 1.9176011588804434e-001; 1.4807083115521585e-001; 9.2364726716986353e-002; 4.5273685465150391e-002; 1.5673473751851151e-002; 3.1554462691875513e-003; 2.3113452403522071e-003; 8.1895392750226735e-004; 2.7524214116785131e-004; 3.5729348198975332e-005; 2.7342206801187888e-006; 2.4676421345798140e-007; 2.1394194479561062e-008; 4.6011760348655917e-010; 3.0972223576062995e-012; 5.4500412650638128e-015; 1.0541326582334014e-018];
  case 19
    n = [0.0000000000000000e+000; 2.4899229757996061e-001; 7.4109534999454085e-001; 1.2304236340273060e+000; 1.7320508075688772e+000; 2.2336260616769419e+000; 2.5960831150492023e+000; 2.8612795760570582e+000; 3.2053337944991944e+000; 3.6353185190372783e+000; 4.1849560176727323e+000; 4.7364330859522967e+000; 5.1870160399136562e+000; 5.6981777684881099e+000; 6.3633944943363696e+000; 7.1221067008046166e+000; 7.9807717985905606e+000; 9.0169397898903032e+000];
    w = [5.1489450806921377e-004; 1.9176011588804437e-001; 1.4807083115521585e-001; 9.2364726716986353e-002; 4.5273685465150523e-002; 1.5673473751851151e-002; 3.1554462691875604e-003; 2.3113452403522050e-003; 8.1895392750226670e-004; 2.7524214116785131e-004; 3.5729348198975447e-005; 2.7342206801187884e-006; 2.4676421345798140e-007; 2.1394194479561056e-008; 4.6011760348656077e-010; 3.0972223576063011e-012; 5.4500412650637663e-015; 1.0541326582337958e-018];
  case 20
    n = [0.0000000000000000e+000; 2.4899229757996061e-001; 7.4109534999454085e-001; 1.2304236340273060e+000; 1.7320508075688772e+000; 2.2336260616769419e+000; 2.5960831150492023e+000; 2.8612795760570582e+000; 3.2053337944991944e+000; 3.6353185190372783e+000; 4.1849560176727323e+000; 4.7364330859522967e+000; 5.1870160399136562e+000; 5.6981777684881099e+000; 6.3633944943363696e+000; 7.1221067008046166e+000; 7.9807717985905606e+000; 9.0169397898903032e+000];
    w = [5.1489450806925551e-004; 1.9176011588804440e-001; 1.4807083115521585e-001; 9.2364726716986298e-002; 4.5273685465150537e-002; 1.5673473751851155e-002; 3.1554462691875573e-003; 2.3113452403522080e-003; 8.1895392750226724e-004; 2.7524214116785137e-004; 3.5729348198975352e-005; 2.7342206801187888e-006; 2.4676421345798124e-007; 2.1394194479561056e-008; 4.6011760348656144e-010; 3.0972223576062963e-012; 5.4500412650638365e-015; 1.0541326582335402e-018];
  case 21
    n = [0.0000000000000000e+000; 2.4899229757996061e-001; 7.4109534999454085e-001; 1.2304236340273060e+000; 1.7320508075688772e+000; 2.2336260616769419e+000; 2.5960831150492023e+000; 2.8612795760570582e+000; 3.2053337944991944e+000; 3.6353185190372783e+000; 4.1849560176727323e+000; 4.7364330859522967e+000; 5.1870160399136562e+000; 5.6981777684881099e+000; 6.3633944943363696e+000; 7.1221067008046166e+000; 7.9807717985905606e+000; 9.0169397898903032e+000];
    w = [5.1489450806913744e-004; 1.9176011588804429e-001; 1.4807083115521594e-001; 9.2364726716986312e-002; 4.5273685465150391e-002; 1.5673473751851151e-002; 3.1554462691875565e-003; 2.3113452403522089e-003; 8.1895392750226670e-004; 2.7524214116785142e-004; 3.5729348198975285e-005; 2.7342206801187888e-006; 2.4676421345798119e-007; 2.1394194479561059e-008; 4.6011760348656594e-010; 3.0972223576062950e-012; 5.4500412650638696e-015; 1.0541326582332041e-018];
  case 22
    n = [0.0000000000000000e+000; 2.4899229757996061e-001; 7.4109534999454085e-001; 1.2304236340273060e+000; 1.7320508075688772e+000; 2.2336260616769419e+000; 2.5960831150492023e+000; 2.8612795760570582e+000; 3.2053337944991944e+000; 3.6353185190372783e+000; 4.1849560176727323e+000; 4.7364330859522967e+000; 5.1870160399136562e+000; 5.6981777684881099e+000; 6.3633944943363696e+000; 7.1221067008046166e+000; 7.9807717985905606e+000; 9.0169397898903032e+000];
    w = [5.1489450806903368e-004; 1.9176011588804448e-001; 1.4807083115521574e-001; 9.2364726716986423e-002; 4.5273685465150516e-002; 1.5673473751851161e-002; 3.1554462691875543e-003; 2.3113452403522063e-003; 8.1895392750226713e-004; 2.7524214116785164e-004; 3.5729348198975319e-005; 2.7342206801187905e-006; 2.4676421345798151e-007; 2.1394194479561082e-008; 4.6011760348656005e-010; 3.0972223576063043e-012; 5.4500412650637592e-015; 1.0541326582339926e-018];
  case 23
    n = [0.0000000000000000e+000; 2.4899229757996061e-001; 7.4109534999454085e-001; 1.2304236340273060e+000; 1.7320508075688772e+000; 2.2336260616769419e+000; 2.5960831150492023e+000; 2.8612795760570582e+000; 3.2053337944991944e+000; 3.6353185190372783e+000; 4.1849560176727323e+000; 4.7364330859522967e+000; 5.1870160399136562e+000; 5.6981777684881099e+000; 6.3633944943363696e+000; 7.1221067008046166e+000; 7.9807717985905606e+000; 9.0169397898903032e+000];
    w = [5.1489450806913755e-004; 1.9176011588804442e-001; 1.4807083115521577e-001; 9.2364726716986381e-002; 4.5273685465150468e-002; 1.5673473751851155e-002; 3.1554462691875560e-003; 2.3113452403522045e-003; 8.1895392750226572e-004; 2.7524214116785158e-004; 3.5729348198975298e-005; 2.7342206801187892e-006; 2.4676421345798129e-007; 2.1394194479561072e-008; 4.6011760348656103e-010; 3.0972223576062963e-012; 5.4500412650638207e-015; 1.0541326582338368e-018];
  case 24
    n = [0.0000000000000000e+000; 2.4899229757996061e-001; 7.4109534999454085e-001; 1.2304236340273060e+000; 1.7320508075688772e+000; 2.2336260616769419e+000; 2.5960831150492023e+000; 2.8612795760570582e+000; 3.2053337944991944e+000; 3.6353185190372783e+000; 4.1849560176727323e+000; 4.7364330859522967e+000; 5.1870160399136562e+000; 5.6981777684881099e+000; 6.3633944943363696e+000; 7.1221067008046166e+000; 7.9807717985905606e+000; 9.0169397898903032e+000];
    w = [5.1489450806914438e-004; 1.9176011588804442e-001; 1.4807083115521577e-001; 9.2364726716986340e-002; 4.5273685465150509e-002; 1.5673473751851155e-002; 3.1554462691875586e-003; 2.3113452403522058e-003; 8.1895392750226551e-004; 2.7524214116785142e-004; 3.5729348198975386e-005; 2.7342206801187884e-006; 2.4676421345798082e-007; 2.1394194479561059e-008; 4.6011760348656382e-010; 3.0972223576062942e-012; 5.4500412650638381e-015; 1.0541326582336941e-018];
  case 25
    n = [0.0000000000000000e+000; 2.4899229757996061e-001; 7.4109534999454085e-001; 1.2304236340273060e+000; 1.7320508075688772e+000; 2.2336260616769419e+000; 2.5960831150492023e+000; 2.8612795760570582e+000; 3.2053337944991944e+000; 3.6353185190372783e+000; 4.1849560176727323e+000; 4.7364330859522967e+000; 5.1870160399136562e+000; 5.6981777684881099e+000; 6.3633944943363696e+000; 7.1221067008046166e+000; 7.9807717985905606e+000; 9.0169397898903032e+000];
    w = [5.1489450806919989e-004; 1.9176011588804437e-001; 1.4807083115521580e-001; 9.2364726716986395e-002; 4.5273685465150426e-002; 1.5673473751851158e-002; 3.1554462691875539e-003; 2.3113452403522054e-003; 8.1895392750226681e-004; 2.7524214116785142e-004; 3.5729348198975292e-005; 2.7342206801187884e-006; 2.4676421345798108e-007; 2.1394194479561056e-008; 4.6011760348655901e-010; 3.0972223576062975e-012; 5.4500412650638412e-015; 1.0541326582337527e-018];
  end

  return
end
