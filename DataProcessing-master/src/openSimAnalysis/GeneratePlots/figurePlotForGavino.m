x = 0:1:100;
fig =figure(2);

VarNames = {'Hip_Flex', 'Hip_Add', 'Hip_Rot' 'Pelvis_Tilt', 'Pelvis_List',...
     'Pelvis_Rot', 'Lumbar_Ext', 'Lumbar_Bend', 'Lumbar_Rot'};

AxisNames = {'Hip Flexion (deg)', 'Hip Adduction (deg)', 'Hip Rotation (deg)', 'Pelvis Tilt (deg)', 'Pelvis List (deg)',...
     'Pelvis Rotation (deg)', 'Trunk Extension (deg)', 'Trunk Bending (deg)', 'Trunk Rotation (deg)'};

% Loop through the structure array and create a subplot for each variable
for n = 1:9
     figure(fig);
     subplot(3,3,n)
     
     % Plot mean waveform data for each of the marker sets separately
     % for the no armour and armour conditions
     h1 = plot(x,importDataCompleteS1.VIRT.(VarNames{n})(:,1)', 'LineWidth',1,'Color','b');
     hold on
     h2 = plot(x, importDataCompleteS1.PHYS.(VarNames{n})(:,1)', 'LineWidth',1,'Color','r');
     hold on
     
     % SD for VIRT marker set
     % Lower bound
     lb = (importDataCompleteS1.VIRT.(VarNames{n})(:,1)...
          - importDataCompleteS1.VIRT.(VarNames{n})(:,2));
     % Upper bound
     ub = (importDataCompleteS1.VIRT.(VarNames{n})(:,1)...
          + importDataCompleteS1.VIRT.(VarNames{n})(:,2));
     % Plot SD
     [ph,msg]=jbfill(x,lb',ub',[0 0 1],[0 0 1],'add',.1);
     
     hold on
     %SD for PHYS marker set
     % Lower bound
     lb = (importDataCompleteS1.PHYS.(VarNames{n})(:,1)...
          - importDataCompleteS1.PHYS.(VarNames{n})(:,2));
     % Upper bound
     ub = (importDataCompleteS1.PHYS.(VarNames{n})(:,1)...
          + importDataCompleteS1.PHYS.(VarNames{n})(:,2));
     [ph,msg]=jbfill(x,lb',ub',[1 0 0],[1 0 0],'add',.1);
     
     box off
     % Place x-labels only for bottom row 
     if n >= 7
          xlabel('Gait Cycle (%)');
     end
     % Include y-axis name for all variables
     ylabel(AxisNames{n})
     hold on
end
% Create legend 
legend([h1 h2], {'New Marker Set', 'Markers on Skin'});