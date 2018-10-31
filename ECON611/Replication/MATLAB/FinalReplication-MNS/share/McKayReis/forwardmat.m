% compute transition matrix
% Pi(i,j) is probability to go from state j to state i!
% modified 9/24/12 by AGM
function Pi = forwardmat(iFullMat,par,R,wage,tau,dividend)
  global Params;

  nd = Params.ndstst;
IR = []; IC = []; VV = [];
  
  
  for ip=1:Params.npp  % this period's state
      
    
    [~,~,Kend] = get_cnbp(Params.knotDistrK,par,R,wage,tau, dividend,ip);

        
    %--
    
    Ti = lineartrans(Params.knotDistrK,Kend);
    ir{ip} = Ti.iTo;  % row indicates to which position we go to
    ic{ip} = Ti.iFr;  % column indicates which position we come from
    vv{ip} = Ti.Val;
    
    
    if(iFullMat)
      offsi = (ip-1)*nd;
      for jp=1:Params.npp  % next period's state

        pp = transProb(ip,jp);
        
        if length(pp) > 1
            pp = reshape(repmat(pp',2,1),2*length(pp),1);
        end
        
        
        if(any(pp>0))
            offsj = (jp-1)*nd;
            IR = [IR;offsj+ir{ip}];  % where to go! take offsi
            IC = [IC;offsi+ic{ip}];  % where come from! take offsi
            VV = [VV;pp.*vv{ip}];  
        end
      end
    else
        Pij{ip} = sparse(ir{ip},ic{ip},vv{ip},nd,nd);
    end
  end  %end loop over this period state
  
  if(iFullMat)
        nn = nd*Params.npp;
        Pi = sparse(IR,IC,VV,nn,nn);
  else
      warning('have not updated this to include variable search effort')
        transMat = Params.Trans
        Pi = op_concat(op_kron(newTransProb',speye(nd)),blockdiag(Pij));  %the ordering here matters for sequence of taking saving decision then learning new skill/employment status rather then vice versa as in Reiter's model.
  end




