function save_markersbox_lake(outpath,k,MTK,MI,MX,MY,MRHO,MFLOW,MMU,MPL,MCP,MKT,MHR,...
            MSXX,MSXY,META,MEXX,MEXY,MPR,MGII,MBII,MRAT,MCXN, MCYN, TXN, plastic_yield, ntimestep)
        
MTK = MTK(k);
MI = MI(k);
MX = MX(k);
MY = MY(k);
MSXX = MSXX(k);
MSXY = MSXY(k);
META = META(k);
MEXX = MEXX(k);
MEXY = MEXY(k);
MPR = MPR(k);
MGII = MGII(k);
MBII = MBII(k);
MRAT = MRAT(k);
MCXN = MCXN(k);
MCYN = MCYN(k);
TXN = TXN(k);
plastic_yield = plastic_yield(k);


save([outpath,'/markers_',num2str(ntimestep),'s.mat'],...
    'MTK','MI','MX','MY','MRHO','MFLOW','MMU','MPL','MCP','MKT','MHR',...
    'MSXX','MSXY','META','MEXX','MEXY','MPR','MGII','MBII','MRAT','MCXN', 'MCYN', 'TXN', 'plastic_yield','ntimestep');