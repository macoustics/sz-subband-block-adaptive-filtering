function FlopCount = GradientDescentComplexity(BufferSize, NumberOfLoudspeakers, NumberOfMicrophones, NumberOfGradientSteps, NumberOfLineSearchSteps, RealSequenceFlag)
% Computational complexity for incomplete gradient descent optimization

% For the new buffer, calculate the spectrum
GradientUpkeep1 = CalculateInputSpectrum(BufferSize, RealSequenceFlag);
% For the new buffer, update the estimate of the input autospectrum
GradientUpkeep2 = UpdateInputAutospectrumEstimate(BufferSize, RealSequenceFlag);
% Total upkeep
GradientUpkeep = GradientUpkeep1 + GradientUpkeep2;

% For every gradient step
% Calculate the gradient
GradientStep1 = EvaluateGradient(BufferSize, NumberOfLoudspeakers, NumberOfMicrophones, RealSequenceFlag);
% Take the gradient step
GradientStep2 = UpdateFilters(BufferSize, NumberOfLoudspeakers, RealSequenceFlag);

% For every gradient step - do a line-search
% For the new line search
% Evaluate the current cost J(w_i)
LineSearchUpkeep1 = EvaluateCostFunction(BufferSize, NumberOfLoudspeakers, NumberOfMicrophones, RealSequenceFlag);
% Evaluate the scaled norm of the gradient
LineSearchUpkeep2 = ScaledGradientNorm(BufferSize, NumberOfLoudspeakers, RealSequenceFlag);
% Total line-search upkeep
LineSearchUpkeep = LineSearchUpkeep1 + LineSearchUpkeep2;

% For every line-search step
% Evaluate potential new filters x = w_i - alpha*Gradient(w_i)
LineSearchStep1 = UpdateFilters(BufferSize, NumberOfLoudspeakers, RealSequenceFlag);
% Evaluate cost at potential new filters J(x)
LineSearchStep2 = EvaluateCostFunction(BufferSize, NumberOfLoudspeakers, NumberOfMicrophones, RealSequenceFlag);
% Evaluate previous scaled cost multiplied by step-size alpha*J(w_i)
LineSearchStep3 = rMult();
% Scale step-size alpha = alpha*rho
LineSearchStep4 = rMult();
% Total for line-search step
LineSearchStep = LineSearchStep1 + LineSearchStep2 + LineSearchStep3 + LineSearchStep4;

% Line-search performed in each gradient step
GradientStep3 = LineSearchUpkeep + NumberOfLineSearchSteps*(LineSearchStep);
% Total gradient step
GradientStep = GradientStep1 + GradientStep2 + GradientStep3;

% Total flop-count for incomplete gradient descent optimization
FlopCount = GradientUpkeep + NumberOfGradientSteps*GradientStep;
end

%% Helper functions
function FlopCount = CalculateInputSpectrum(BufferSize, RealSequenceFlag)
% Calculate U[s]
% Exploit half-spectrum symmetry if possible
if RealSequenceFlag
    BufferSize = BufferSize/2;
end

FlopCount = FourierTransform(BufferSize, RealSequenceFlag);
end

function FlopCount = UpdateInputAutospectrumEstimate(BufferSize, RealSequenceFlag)
% P[s] = y*P[s-1] + (1-y)*U^H[s] * U[s]
% P[s] = diag(p[s]) \in \mathbb C^{BufferSize x BufferSize}
% U[s] = diag(u[s]) \in \mathbb C^{BufferSize x BufferSize}
% Assume u[s], P[s-1], and (1-y) to be precomputed

% Exploit half-spectrum symmetry if possible
if RealSequenceFlag
    BufferSize = BufferSize/2+1;
end

% x1 = U^H[s] * U[s] = diag(abs(u[s]).^2)
CostX1 = BufferSize*cAbsSq();
% x2 = (1-y)*x1
CostX2 = BufferSize*rMult();
% x3 = y*P[s-1]
CostX3 = BufferSize*rMult();
% x4 = x3 + x2
CostX4 = BufferSize*rAdd();

FlopCount = CostX1 + CostX2 + CostX3 + CostX4;
end

function FlopCount = EvaluateSingleMicrophoneError(BufferSize, NumberOfLoudspeakers, RealSequenceFlag)
% e^{m}[s] = DFT Z^{01} IDFT U[s](R^{m,t}w^{t} - R^{m} w)
% R^{m} = [R^{m,1}, R^{m,2}, ..., R^{m,NumberOfLoudspeakers}]
% R^{m,l} = diag(r^{m,l}) \in \mathbb C^{BufferSize x BufferSize}
% R^{m,t} = diag(r^{m,t}) \in \mathbb C^{BufferSize x BufferSize}
% w^{t} \in \mathbb C^{BufferSize}
% w \in \mathbb C^{BufferSize*NumberOfLoudspeakers}
% U[s] = diag(u[s]) \in \mathbb C^{BufferSize x BufferSize}
% Assume u[s], R^{m}, R^{m,t}, and w^{t} to be precomputed

% Exploit half-spectrum symmetry if possible
if RealSequenceFlag
    BufferSize = BufferSize/2+1;
end

% x1 = R^{m} w
CostX1 = BufferSize*NumberOfLoudspeakers*cMult();
% x2 = R^{m,t} w^{t}
CostX2 = BufferSize*cMult();
% CostX2 = BufferSize*NumberOfLoudspeakers*cMult();
% x3 = x2 - x1
CostX3 = BufferSize*cAdd();
% x4 = U[s] * x3
CostX4 = BufferSize*cMult();
% x5 = IDFT x4
CostX5 = FourierTransform(BufferSize, RealSequenceFlag);
% x6 = Z^{01} x5 (set leading BufferSize/2 values to 0)
CostX6 = 0;
% x7 = DFT x6
CostX6 = FourierTransform(BufferSize, RealSequenceFlag);

FlopCount = CostX1 + CostX2 + CostX3 + CostX4 + CostX5 + CostX6;
end

function FlopCount = EvaluateGradient(BufferSize, NumberOfLoudspeakers, NumberOfMicrophones, RealSequenceFlag)
% grad = kron(I_L, DFT Z^{10} IDFT) (lambda*kron(I_L, P[s]) w -
% sum_mics(R^{m H} U^H[s] e^m[s])
% e^{m}[s] \in \mathbb C^{BufferSize}
% R^{m} = [R^{m,1}, R^{m,2}, ..., R^{m,NumberOfLoudspeakers}]
% R^{m,l} = diag(r^{m,l}) \in \mathbb C^{BufferSize x BufferSize}
% w \in \mathbb C^{BufferSize*NumberOfLoudspeakers}
% U[s] = diag(u[s]) \in \mathbb C^{BufferSize x BufferSize}
% I_L = eye(NumberOfLoudspeakers)
% Assume u[s], R^{m}, P[s], and w^{t} to be precomputed

% Exploit half-spectrum symmetry if possible
if RealSequenceFlag
    BufferSize = BufferSize/2+1;
end

% x1 = lambda*P[s]
CostX1 = BufferSize*rMult();
% x2 = kron(I_L,x1) w
CostX2 = BufferSize*NumberOfLoudspeakers*rcMult();
% x3a^{m} = EvaluateSingleMicrophoneError
CostX3a = EvaluateSingleMicrophoneError(BufferSize, NumberOfLoudspeakers, RealSequenceFlag);
% x3b^{m} = U^H[s] x3a^{m}
CostX3b = BufferSize*cMult();
% x3c^{m} = R^{m H} x3b^{m}
CostX3c = BufferSize*NumberOfLoudspeakers*cMult();
% x3d = sum_mics(x3c^{m})
CostX3d = (NumberOfMicrophones-1)*cAdd();
% Combined cost for x3
CostX3 = NumberOfMicrophones*(CostX3a + CostX3b + CostX3c) + CostX3d;
% x4 = x2 - x3
CostX4 = BufferSize*NumberOfLoudspeakers*cAdd();
% x5 = kron(I_L, DFT Z^{10} IDFT) x4
% x5a = IDFT x4(1:BufferSize)
CostX5a = FourierTransform(BufferSize, RealSequenceFlag);
% x5b = Z^{10} x5a (set lagging BufferSize/2 values to 0)
CostX5b = 0;
% x5c = DFT x5b
CostX5c = FourierTransform(BufferSize, RealSequenceFlag);
% Combined cost for x5
CostX5 = NumberOfLoudspeakers*(CostX5a + CostX5b + CostX5c);

FlopCount = CostX1 + CostX2 + CostX3 + CostX4 + CostX5;
end

function FlopCount = EvaluateCostFunction(BufferSize, NumberOfLoudspeakers, NumberOfMicrophones, RealSequenceFlag)
% cost = lambda*w^H kron(I_L, P[s]) w + sum_mics(e^{m H}[s] e^m[s])
% e^{m}[s] \in \mathbb C^{BufferSize}
% w \in \mathbb C^{BufferSize*NumberOfLoudspeakers}
% P[s] = diag(p[s]) \in \mathbb C^{BufferSize x BufferSize}
% I_L = eye(NumberOfLoudspeakers)
% Assume P[s] and w to be precomputed

% Exploit half-spectrum symmetry if possible
if RealSequenceFlag
    BufferSize = BufferSize/2+1;
end

% x1 = abs(w).^2
CostX1 = BufferSize*NumberOfLoudspeakers*cAbsSq();
% x2 = kron(I_L,P[s]) x1
CostX2 = BufferSize*NumberOfLoudspeakers*rMult();
% x3 = sum(x2)
CostX3 = (BufferSize*NumberOfLoudspeakers-1)*rAdd();
% x4 = lambda*x3
CostX4 = rAdd();
% x5a^{m} = EvaluateSingleMicrophoneError
CostX5a = EvaluateSingleMicrophoneError(BufferSize, NumberOfLoudspeakers, RealSequenceFlag);
% x5b^{m} = abs(x5a).^2
CostX5b = BufferSize*cAbsSq();
% x5c^{m} = sum(x5b^{m})
CostX5c = (BufferSize-1)*rAdd();
% x5d = sum_mics(x5c^{m})
CostX5d = (NumberOfMicrophones-1)*rAdd();
% Combined cost for x5
CostX5 = NumberOfMicrophones*(CostX5a + CostX5b + CostX5c) + CostX5d;
% cost = x6 = x4 + x5
CostX6 = rAdd();

FlopCount = CostX1 + CostX2 + CostX3 + CostX4 + CostX5 + CostX6;
end

function FlopCount = UpdateFilters(BufferSize, NumberOfLoudspeakers, RealSequenceFlag)
% w_{i+1} = w_{i} - alpha*Gradient
% w \in \mathbb C^{BufferSize*NumberOfLoudspeakers}
% Assume Gradient and w to be precomputed

% Exploit half-spectrum symmetry if possible
if RealSequenceFlag
    BufferSize = BufferSize/2+1;
end

% x1 = -alpha*Gradient
CostX1 = BufferSize*NumberOfLoudspeakers*rcMult();
% x2 = w_{i} + x1
CostX2 = BufferSize*NumberOfLoudspeakers*cAdd();

FlopCount = CostX1 + CostX2;
end

function FlopCount = ScaledGradientNorm(BufferSize, NumberOfLoudspeakers, RealSequenceFlag)
% c_1 * norm(Gradient)
% Gradient \in \mathbb C^{BufferSize*NumberOfLoudspeakers}
% Assume Gradient to be precomputed

% Exploit half-spectrum symmetry if possible
if RealSequenceFlag
    BufferSize = BufferSize/2+1;
end

% x1 = abs(Gradient).^2
CostX1 = BufferSize*NumberOfLoudspeakers*cAbsSq();
% x2 = sum(x1)
CostX2 = (BufferSize*NumberOfLoudspeakers-1)*rAdd();
% x3 = c_1*x2
CostX3 = rMult();

FlopCount = CostX1 + CostX2 + CostX3;
end