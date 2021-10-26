function optimNames = LoadOptimNames()
    % FMC <Common struct optimNames because the number of names is not very high>
    optimNames.OptimizationFMCName = 'fmcconfig';
    optimNames.OptimizationFMCUseName = 'use_fmc';
    optimNames.OptimizationFMCUStartName = 'start_from';
    optimNames.OptimizationFMCX0Name = 'x0';
    optimNames.OptimizationFMCRestartsName = 'num_starts';
    % TLBO
    optimNames.OptimizationTLBOName = 'tlboconfig';
    optimNames.OptimizationTLBOUseName = 'use_tlbo';
    optimNames.OptimizationTLBOPopSizeName = 'pop_size';
    optimNames.OptimizationTLBONCyclesName = 'nun_cycles';
    optimNames.OptimizationTLBORestartsName = 'num_starts';
    % MUMSA
    optimNames.OptimizationMUMSAName = 'mumsaconfig';
    optimNames.OptimizationMUMSAUseName = 'use_mumsa';
    optimNames.OptimizationMUMSANPName = 'np';
    optimNames.OptimizationMUMSACPName = 'cp';
    optimNames.OptimizationMUMSAMPName = 'mp';
    optimNames.OptimizationMUMSARangeName = 'range';
    optimNames.OptimizationMUMSAFName = 'f';
    optimNames.OptimizationMUMSAItermaxName = 'itermax';
    optimNames.OptimizationMUMSARestartsName = 'num_starts';
    % DE
    optimNames.OptimizationDEName = 'deconfig';
    optimNames.OptimizationDEUseName = 'use_de';
    optimNames.OptimizationDENPName = 'np';
    optimNames.OptimizationDECRName = 'cr';
    optimNames.OptimizationDEFName = 'f';
    optimNames.OptimizationDEXName = 'x';
    optimNames.OptimizationDEYName = 'y';
    optimNames.OptimizationDEGenMaxName = 'gen_max';
    optimNames.OptimizationDERestartsName = 'num_starts';
    
    %----------RESULT NAMES--------
    optimNames.fmcResultsName = 'fmcResults';
    optimNames.tlboResultsName = 'tlboResults';
    optimNames.mumsaResultsName = 'mumsaResults';
    optimNames.deResultsName = 'deResults';
    
    optimNames.bestRawResultName = 'bestRawResult';
    optimNames.bestRawFuncValName = 'bestRawFuncVal';
    optimNames.meanQualityName = 'meanQuality';
    optimNames.stdQualityName = 'stdQuality';
    optimNames.meanTimeName = 'meanTime';
    
    optimNames.bestXName = 'X1';
    optimNames.bestXRawValName = 'X1RawVal';
    optimNames.UpdatedTh2Name = 'th2_prescrib';
    optimNames.XYachievedName = 'XYachieved';
    optimNames.costAEName = 'costAE';
    optimNames.costNSDVName = 'costNSDV';
    optimNames.costTDName = 'costTD';
    optimNames.LogName = 'log';
end
