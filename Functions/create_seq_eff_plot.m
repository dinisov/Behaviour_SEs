%this plotting function will hopefully make all other seq eff plotting functions redundant
function create_seq_eff_plot(RT_data,fit_data,varargin)

%parse optional arguments
names = {'xlim' 'ylim' 'font' 'labeltype' 'pdfready' 'selabels' 'markersize' 'linewidth' 'fontsize' 'histlength' 'errors' 'scores','reOrder','n_back'};%options
dflts = {[] [] 'times' 'ar' 'yes' 'on' 20 2 16 4 [] [],[1:16],5};%defaults
[xl,yl,font,labeltype,pdfready,selabels,markersize,linewidth,fontsize,histlength,errors,scores, reOrder,n_back] = parseArgsUser(names, dflts, 'userargs', varargin{:});

hold all;

%Check if reordering requested
%if ~isequal(reOrder,[1:16])
%    disp(['-- SEs reordering requested --'])
%end
if numel(unique(reOrder)) ~= size(RT_data,1)
    ['## ALERT: REQUESTED REORDERING OMITS OR OVERLENGTHS RECEIVED DATA ##']
    ['Proposed order: ',num2str(reOrder)]
    ['Received # of sequences: ',num2str(size(RT_data,1))]
end
if 2^histlength ~= size(RT_data,1)
    ['## Alert: Disparity between histlength (',num2str(histlength),' [-> ',num2str(2^histlength),' seqs]) and sequences in data [',num2str(size(RT_data,1)),'] ##']
    crash = yes
end

if ~isempty(errors)
    
    if size(RT_data,2) == size(errors,2) %if simple SEM or SE are to be displayed
    
        for i = 1:size(RT_data,2)
            %errorbar(1:2^histlength,RT_data(:,i),errors(:,i),'marker','.','markersize',markersize,'linewidth',linewidth);
            errorbar(1:2^histlength,RT_data([reOrder],i),errors([reOrder],i),'marker','.','markersize',markersize,'linewidth',linewidth); %Note: Will have extreme unintended effects if proposed reordering not same size as actual data
        end
        xlim([0,2^histlength+1])
        xticks([1:2^histlength])
    
        data_min = min(min([RT_data-errors fit_data])); 
        data_max = max(max([RT_data+errors fit_data]));
        
    else%if there are different upper and lower bounds CIs

        %Check if reordering requested
        if ~isequal(reOrder,[1:16])
            disp(['-# Caution: SEs reorder requested but inapplicable to custom error; Not applying #-'])
        end
        
        for i = 1:size(RT_data,2)
            errorbar(1:2^histlength,RT_data(:,i),errors(:,2*i-1),errors(:,2*i),'marker','.','markersize',markersize,'linewidth',linewidth);
        end
        
        data_min = min(min([RT_data-errors(:,1:2:size(errors,2)) fit_data])); 
        data_max = max(max([RT_data+errors(:,2:2:size(errors,2))  fit_data]));
        
    end
    
else

    for i = 1:size(RT_data,2)
        %plot(RT_data(:,i),'marker','.','markersize',markersize,'linewidth',linewidth); 
        plot(RT_data([reOrder],i),'marker','.','markersize',markersize,'linewidth',linewidth); 
    end

    data_min = min(min([RT_data fit_data])); data_max = max(max([RT_data fit_data]));
    
end

data_range = data_max-data_min;

if exist('fit_data','var')
    hold on; plot(fit_data,'r--','marker','+','markersize',12,'linewidth',linewidth);
end
        
if data_range ~= 0

    if isempty(xl)
        xlim([.5 2^histlength+.5]);
    else
        xlim(xl);
    end

    if isempty(yl) && data_range ~= 0
            ylim([data_min-0.02*data_range data_max+0.02*data_range]);       
    else
            ylim(yl);  
    end

end

set(gca,'xtick',[],'Units','Pixels');

%make x labels as separate text strings
if strcmpi('on',selabels)
    
    add_se_labels(font,labeltype,fontsize,histlength, reOrder, n_back);
    
end

%remove all possible white from picture; wrap figure around labels (not automatic since they are text strings)
%ready to print via printpdf()
if strcmpi('yes',pdfready) && strcmpi('on',selabels)
    if n_back <= 5
        axis square
    end

    set(gca,'units','centimeters'); set(gcf,'units','centimeters');  set(gca,'fontsize',14);

    %absolute positions, not very elegant but necessary in the case of multiple lines
    %there is no immediate way to translate between the units of the text box (data) and axes/figure units
    switch labeltype
        case 'ar'
            if n_back <= 5
                set(gcf,'position',[11.5 7.25 12 12.4]); %Absolute X, Y, W, H
                set(gca,'Position',[1.35 2.73 9.75 8.28]);
            elseif n_back == 6
                set(gcf,'position',[11.5 5.25 16.25 15]); %Absolute X, Y, W, H
                set(gca,'Position',[1.35 3.73 12.75 8.28]);
            else
                set(gcf,'units','normalized','outerposition',[0 0 1 1]) %Fullscreen figure itself
                set(gca,'units','normalized','outerposition',[0 0+0.1*(n_back-5) 1 1-0.1*(n_back-5)]) %Empirical calcs to show all legend/etc
            end
        case 'xy'
            set(gcf,'position',[11.76 6.83 10.1 11.7]);
            set(gca,'Position',[1.04 3.3 9.74 8.28]);
         case 'lr'
            set(gcf,'position',[11.5 7.25 10.4 11.93]);
            set(gca,'Position',[1.44 2.83 9.75 8.28]);
         case 'ci'
            set(gcf,'position',[11.5 7.25 10.1 11.28]);
            set(gca,'Position',[1.1 2.66 9.8 8.28]);
    end
    
end

if ~isempty(scores)
   
    axes('Units','normalized','Position',[.8 .8 .15 .15]);
    h = bar(1:3,diag(scores),'stacked');
    set(h(1),'facecolor','red'); set(h(2),'facecolor','green'); set(h(3),'facecolor','blue')
    
end

end

function add_se_labels(font,labeltype,fontsize,histlength, reOrder, n_back)
    
    currentYlim = get(gca,'ylim'); ylim_range = range(currentYlim);
    %x_labels = choose_labels(font,labeltype, reOrder);   
    x_labels = choose_labels(font,labeltype, reOrder, n_back);
    
    %this automatically adjusts for any length of A/R history from 1 to 4 (2 to 5 in terms of X/Y)
    for i = 1:2^histlength
        %text(i,currentYlim(1)-.01*ylim_range,x_labels{i*2^(-histlength+4)}((-histlength+5):end,:),'VerticalAlignment','top','HorizontalAlignment',...
        %'center','fontsize',fontsize,'interpreter','latex'); %Hardcoded 5-back
        text(i,currentYlim(1)-.01*ylim_range,x_labels{i*2^(-histlength+(n_back-1))}((-histlength+n_back):end,:),'VerticalAlignment','top','HorizontalAlignment',...
        'center','fontsize',fontsize,'interpreter','latex');
    end
    
end

function [labels] = choose_labels(font,labeltype, reOrder, n_back)

    switch n_back
        case 5
    
             switch font   
                case 'times'
                    switch labeltype
                        case 'ar'
                            load binomial_x_labels_latex_alt_rep
                        case 'xy'
                            load binomial_x_labels_latex_x_y
                        case 'lr'
                            load x_labels_latex_left_right
                        case 'ci'
                            load x_labels_latex_coherent_incoherent                  
                    end
                case 'sans'   
                    switch labeltype
                        case 'ar'
                            load binomial_x_labels_latex_alt_rep_sans
                        case 'xy'
                            load binomial_x_labels_latex_x_y_sans
                        case 'lr'
                            load x_labels_latex_left_right_sans
                        case 'ci'
                            load x_labels_latex_coherent_incoherent_sans                          
                    end
             end
             
             labels = binomial_x_labels_latex;

        otherwise
            loadName = [num2str(n_back),'-back_legend.mat'];
            %load 6-back_legend.mat
            eval(['load ',loadName]) %Return of the evil eval

            %No font specification

            switch labeltype
                case 'ar'
                    %labels = anynomial_x_labels_latex;
                    labels = anynomial_x_labels_latex_canonical;
                case 'xy'
                    %labels = anynomial_legend;
                    labels = anynomial_legend_canonical;
                case 'orig' %New
                    %labels = anynomial_legend_orig;
                    labels = anynomial_legend_orig_canonical;
            end

    end

     %QA
     if numel(unique(reOrder)) ~= size(labels,2)
         ['## Alert: Requested reordering differs in apparent size from supplied labels ##']
         crash = yes
     end
     labels = labels(reOrder); %Under default circumstances won't change anything
     %if ~isequal(reOrder,[1:16]) %Hardcoded default
     %    disp(['-# Labels reordered #-'])
     %end
     
end