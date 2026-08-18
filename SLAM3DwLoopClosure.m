rng(0);
bag=rosbag('C:\Users\ACER\Downloads\HolybroStdn01.bag'); %update filepath according to your device
lidarBag = select(bag, 'Topic', '/ouster/points');
nFrames=lidarBag.NumMessages;
pGraph = poseGraph3D;
[maxLdrRng, refVector, maxDist, maxAngDist, gridStep, distMoveThres] = deal(20,[0,0,1],0.15,15,1.5,0.3);
[loopCloseSearchRadius, subMapThres,rmseThres,optIntr]=deal(8,50,0.45,2);
infoMat = [1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,1];
nLoopClosures = 0; 
mapUpt = false;
[prevPc,prevTform]=deal([],[]);
[scanAccepted,j] = deal(0,0);
omap = occupancyMap3D(15);
pcProcessed = cell(1,nFrames);
pcsToView = cell(1,nFrames); 
for i=1:nFrames
    msg = readMessages(lidarBag, i, 'DataFormat', 'struct');
    curPtCloud=pointCloud(rosReadXYZ(msg{1}));
    xyzMatrix=curPtCloud.Location ;

    ind=(-maxLdrRng<xyzMatrix(:,1) & maxLdrRng>xyzMatrix(:,1) & -maxLdrRng<xyzMatrix(:,2) & maxLdrRng>xyzMatrix(:,2) & (abs(xyzMatrix(:,2)) > abs(1.2 * xyzMatrix(:,1)) | xyzMatrix(:,1) > 0));
    pcl=pointCloud(xyzMatrix(ind,:));

    [~,~,outliers]=pcfitplane(pcl,maxDist,refVector,maxAngDist);
    pclXgnd=select(pcl,outliers,'OutputSize','full');
    [~,~,outliers]=pcfitplane(pclXgnd,0.2,refVector,maxAngDist);
    pclXgnd=select(pclXgnd,outliers,'OutputSize','full');
    ind=((pclXgnd.Location(:,3)<0.4185)&(pclXgnd.Location(:,3)>-0.6871));
    pclFiltered = select(pclXgnd, ind, 'OutputSize', 'full');
    pclFinl = pcdownsample(pclFiltered,'random',.25);

    if j==0
        tform=[];
        scanAccepted=1;
    else
        if j==1
            tform=pcregisterndt(pclFinl,prevPc,gridStep);
        else
            tform=pcregisterndt(pclFinl,prevPc,gridStep,'InitialTransform',prevTform);
        end
        relPose=[tform2trvec(tform.T') tform2quat(tform.T')];
        if norm(relPose(1:3))>distMoveThres
            addRelativePose(pGraph,relPose);
            scanAccepted=1;
        else
            scanAccepted=0;
        end
    end
    if scanAccepted == 1
        j=j+1;
        pcProcessed{j}=pclFinl;
        if j>subMapThres
            nodesall=nodes(pGraph);
            curxyz=nodesall(j,1:3);
            prevposmat=nodesall(1:j-subMapThres,1:3);
            dists=vecnorm(prevposmat-curxyz,2,2);
            probcand=find(dists<loopCloseSearchRadius);
        
            rmseMin=inf;
            for k=1:length(probcand)
                [lptform,~,rmse]=pcregisterndt(pclFinl,pcProcessed{probcand(k)},gridStep);
                if rmse<rmseMin
                    rmseMin=rmse;
                    bstTform=lptform;
                    loopcand=probcand(k);
                end
            end
            if rmseMin<rmseThres
                relPose=[tform2trvec(bstTform.T') tform2quat(bstTform.T')];
                addRelativePose(pGraph,relPose,infoMat,loopcand,j);
                nLoopClosures=nLoopClosures+1;
            end
        end
        if (nLoopClosures==optIntr)||((nLoopClosures>0)&&(i==nFrames))
            if loopCloseSearchRadius~=1 
                disp('doing PGO to reduce drift'); 
            end
            pGraph=optimizePoseGraph(pGraph);
            nLoopClosures=0;
            optIntr=optIntr*7;
        end
        pcsToView{j}=pclFinl;
        if rem(j,15)==0 
            fprintf('.'); end
    


    end
    if ~isempty(pclFinl) 
        prevPc=pclFinl; end
    if scanAccepted==1 && ~isempty(tform) 
        prevTform=tform; end

    clear msg curPtCloud xyzMatrix pcl pclXgnd ind pclFiltered lptform;
end

            
 
% --- Final Output Generation (Native Rendering) ---
disp('Building Final 3D Map...');
nodesPos=nodes(pGraph);

for k=1:(size(nodesPos,1)-1)
    insertPointCloud(omap,nodesPos(k,:),pcsToView{k}.removeInvalidPoints,maxLdrRng);
end

% Native MATLAB Rendering (Replaces the Helper Function)
figure('Name','Final trajectory','Color','w');
show(omap); 
hold on;
show(pGraph,'IDs','off'); 
title('Optimized 3D Occupancy Map');
grid on; axis equal; view(3);     