function [A,InvNum,InvList,aviobj]=InvPFuncListed(A,StartPoint,BW,DensDef,DensInv,Grav,CapPres,vidnovid)

InvNum = max(A(:))+1; %Value to set at invaded pore space
i = StartPoint(1);
j = StartPoint(2);
A(i,j) = InvNum;
it = 1; %iteration counter

if vidnovid == 1
    aviobj=avifile('flowingTEMP.avi');
    aviobj.fps=100;
    aviobj.quality=50;
    aviobj.compression='Cinepak';
end

while BW < i < size(A,1) && BW < j < size(A,2)
    B = [A(i-1,j-1) A(i-1,j) A(i-1,j+1); A(i,j-1) A(i,j) A(i,j+1); A(i+1,j-1) A(i+1,j) A(i+1,j+1)]; 
    if B(1,2) ~= InvNum
        B(1,2) = B(1,2) + (DensDef-DensInv)*-Grav*((size(A,2) - (j-1))*1e-3); %if location is not invaded, compute Ptotal
    end
    if B(1,1) ~= InvNum
        B(1,1) = B(1,1) + 0.5*(DensDef-DensInv)*-Grav*((size(A,2) - (j-1))*1e-3);
    end
    if B(1,3) ~= InvNum
        B(1,3) = B(1,3) + 0.5*(DensDef-DensInv)*-Grav*((size(A,2) - (j-1))*1e-3);
    end
    if B(3,2) ~= InvNum
        B(3,2) = B(3,2) + (DensDef-DensInv)*Grav*((size(A,2) - (j+1))*1e-3); %note gravity sign change
    end
    if B(3,1) ~= InvNum
        B(3,1) = B(3,1) + 0.5*(DensDef-DensInv)*Grav*((size(A,2) - (j+1))*1e-3);
    end
    if B(3,3) ~= InvNum
        B(3,3) = B(3,3) + 0.5*(DensDef-DensInv)*Grav*((size(A,2) - (j+1))*1e-3);
    end
    
    InvList((8*it-7):(8*it),1) = [i-1 i-1 i-1 i i i+1 i+1 i+1];
    InvList((8*it-7):(8*it),2) = [j-1 j  j+1 j-1 j+1 j-1 j j+1];
    InvList((8*it-7):(8*it),3) = [B(1,1) B(1,2) B(1,3) B(2,1) B(2,3) B(3,1) B(3,2) B(3,3)];
    
    mlist = min(InvList(:,3)); 
    mrow = find(InvList(:,3)==mlist);
    mloc = InvList(mrow(1),1:2); %mrow(1) because mrow will find repeated entries in mlist
    i = mloc(1);
    j = mloc(2);
    if i < BW || j < BW || i > size(A,1)-BW || j > size(A,2)-BW
        break
    end
    A(i,j) = InvNum;
    InvList(mrow,3)= InvNum;
    it = it+1;
    
    if vidnovid == 1
        PathSub = A - CapPres;
        Path = PathSub>=InvNum;
        imshow(Path)
        F=getframe;
        aviobj=addframe(aviobj,F);
    end
end

if vidnovid == 1 
    aviobj=close(aviobj);
end
