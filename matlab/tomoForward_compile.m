% TOMOFORWARD_SCRIPT   Generate dynamic library tomoForward from tomoForward.
% 
% Script generated from project 'tomoForward.prj' on 27-Mar-2024.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.EmbeddedCodeConfig'.
cfg = coder.config('dll','ecoder',true);
cfg.HardwareImplementation.ProdHWDeviceType = 'Google->V8 Engine';
cfg.HardwareImplementation.TargetHWDeviceType = 'Google->V8 Engine';
cfg.TargetLang = 'C++';
cfg.GenerateReport = true;
cfg.MaxIdLength = 1024;
cfg.PostCodeGenCommand = 'wasm.coder.postcodegen.addCLinkage(buildInfo); wasm.coder.postcodegen.registerExportedFunctions(buildInfo)';
cfg.ReportPotentialDifferences = false;
cfg.Toolchain = 'Emscripten v3.1.56 | gmake (64-bit Windows)';

%% Define argument types for entry-point 'tomoForward'.
ARGS = cell(1,1);
ARGS{1} = cell(4,1);
ARGS{1}{1} = coder.typeof(0);
ARGS{1}{2} = coder.typeof(0);
ARGS{1}{3} = coder.typeof(0);
ARGS{1}{4} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg tomoForward -args ARGS{1}

