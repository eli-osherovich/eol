function test_suite = testUnitaryDFT
initTestSuite;


% Test that the constructor accepts only one or two inputs
function testNArgIn
% The constructor requires either one or two input artuments (positive
% intergers).

% Try zero inputs.
f = @() UnitaryDFT_eo();
assertExceptionThrown(f, 'EOL:UnitaryDFT:WrongArgNum');

% Try 3 inputs.
f = @() UnitaryDFT_eo(1, 2, 3);
assertExceptionThrown(f, 'MATLAB:maxrhs');


function testForward
% Without zero-padding.
in = rand(2,3) + 1i*rand(2,3);
expectedOut = fftn(in)/sqrt(numel(in));
A = UnitaryDFT_eo(size(in));
assertElementsAlmostEqual(expectedOut(:), A*in);

% With zero-padding
in = rand(2,3) + 1i*rand(2,3);
expectedOut = fftn(in, [4, 6])/sqrt(4*6);
A = UnitaryDFT_eo(size(in), [4, 6]);
assertElementsAlmostEqual(expectedOut(:), A*in);

function testAdjoint
% Without zero-padding.
in1 = rand(2,3) + 1i*rand(2,3);
in2 = rand(2,3) + 1i*rand(2,3);
A = UnitaryDFT_eo(size(in1));
res1 = in2(:)'*(A*in1);
res2 = (A'*in2)'*in1(:);
% Verify that the operator is indeed adjoint.
assertElementsAlmostEqual(res1, res2);


% With zero-padding.
in1 = rand(2,3) + 1i*rand(2,3);
in2 = rand(4,8) + 1i*rand(4,8);
A = UnitaryDFT_eo(size(in1), size(in2));
res1 = in2(:)'*(A*in1);
res2 = (A'*in2)'*in1(:);
% Verify that the operator is indeed adjoint.
assertElementsAlmostEqual(res1, res2);

