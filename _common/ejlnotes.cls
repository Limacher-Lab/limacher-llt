% aeronautical course notes class
% Save as: aeronotes.cls
% Engine: pdfLaTeX
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{ejlnotes}[2025/08/08 v1.0 Course notes class]

\typeout{=== USING EJLNOTES FROM: ejlnotes.cls ===}

% ----- Base class -----
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax
\LoadClass[11pt,letter]{article}

% ----- Packages (pdfLaTeX friendly) -----
\RequirePackage[T1]{fontenc}
\RequirePackage[utf8]{inputenc}
\RequirePackage{mathptmx}                 % Times-like serif body and math
\RequirePackage[scaled=0.92]{helvet}      % Helvetica-like sans for headings
\RequirePackage{courier}                  % Monospaced text
\RequirePackage{sectsty}
\RequirePackage{titlesec}
\RequirePackage{geometry}
\RequirePackage{fancyhdr}
\RequirePackage{xcolor}
\RequirePackage{graphicx}
\RequirePackage{amsmath,amssymb}
\RequirePackage{microtype}
\RequirePackage[hidelinks]{hyperref}
\RequirePackage{lastpage}
\RequirePackage[nameinlink,capitalise,noabbrev]{cleveref}
\RequirePackage[most]{tcolorbox}

% ----- Geometry (users may override in document) -----
\geometry{left=25mm,right=25mm,top=30mm,bottom=30mm}

% ----- Defaults & user-configurable fields -----
\makeatletter
%\newcommand{\coursename}{\@latex@warning@no@line{\string\coursename not set}}
%\newcommand{\lecttopic}{\@latex@warning@no@line{\string\lecttopic not set}}
%\newcommand{\authorname}{\@latex@warning@no@line{\string\authorname not set}}
\newcommand{\coursename}{}
\newcommand{\lecttopic}{}
\newcommand{\authorname}{}
\makeatother
\newcommand{\datemodified}{\today}

% Font aliases (user may \renewcommand in document)
\newcommand{\headerfont}{\sffamily\bfseries}
\newcommand{\titlefont}{\Large\bfseries}
\newcommand{\bodyfont}{\normalfont}

% math formatting
\RequirePackage{bm}         % bold math, incl. Greek
\newcommand{\vecsym}[1]{\bm{#1}}                 % vectors: bold italic
\newcommand{\tensym}[1]{\bm{\mathsf{#1}}}        % tensors: bold sans-serif

% ----- Section styling -----
\allsectionsfont{\sffamily\bfseries}
\titlespacing*{\section}{0pt}{1ex}{0ex}

%% For date and time
\let\latexToday\today
\RequirePackage[]{datetime}
\let\today\latexToday

% ----- Header & footer -----
\setlength{\headheight}{18pt}
\setlength{\headsep}{10mm}
\pagestyle{fancy}
\fancyhf{}
%\fancyhead[L]{\headerfont \small \coursename\\\headerfont\Large \lecttopic}
%\fancyhead[R]{\headerfont \small By: \authorname\\\headerfont\small Last modified: \datemodified}
\fancyhead[L]{\headerfont\shortstack[l]{\small \coursename\\\Large \lecttopic}}
\fancyhead[R]{\headerfont\shortstack[r]{\small By: \authorname\\\small Last modified: \datemodified}}
\renewcommand{\headrulewidth}{0.8pt}
\renewcommand{\headrule}{\hbox to\headwidth{\leaders\hrule height \headrulewidth\hfill}}
\fancyfoot[C]{\thepage\ of \pageref{LastPage}}

% ----- Paragraph formatting -----
\setlength{\parindent}{0pt}
\setlength{\parskip}{6pt}

% ---- List formatting ------
\RequirePackage{enumitem}
\setlist{
  itemsep=0.3ex,
  topsep=0.8ex,
  parsep=2pt,
  partopsep=0pt,
  leftmargin=5ex
}

% ----- Callout boxes: numbered & referenceable -----
\tcbuselibrary{theorems}
% Counter reset by section, printed as <section>.<n>
\newcounter{callout}[section]
%\renewcommand{\thecallout}{\thesection.\arabic{callout}}
\renewcommand{\thecallout}{\arabic{callout}}
\crefname{callout}{box}{boxes}
\Crefname{callout}{Box}{Boxes}

% Base style (users can override via \tcbset in document)
\tcbset{calloutstyle/.style={
		enhanced, breakable,
		colback=black!2, colframe=black!25,
		boxrule=0.4pt, arc=2pt,
		left=1em, right=1em, top=0.6em, bottom=0.6em,
		borderline west={2pt}{0pt}{black!50},
		fonttitle=\bfseries\sffamily,
		before upper={\setlength{\parskip}{0.6\baselineskip}\setlength{\parindent}{0pt}}
}}

% Untitled callout: "Box <num>"
\newenvironment{callout}[1][]% optional tcolorbox keys
{\refstepcounter{callout}%
	\begin{tcolorbox}[calloutstyle,title={Box~\thecallout},#1]}%
	{\end{tcolorbox}}

% Titled callout: "Box <num>  <title>"
\newenvironment{callouttitle}[2][]% [tcb-keys]{Title}
{\refstepcounter{callout}%
	\begin{tcolorbox}[calloutstyle,title={Box~\thecallout\quad #2},#1]}%
	{\end{tcolorbox}}

% Floating variants (wrap in figure floats)
\RequirePackage{float}
\newenvironment{calloutfloat}[2][tbp]% [placement]{Title}
{\refstepcounter{callout}%
	\begin{figure}[#1]\centering
		\begin{tcolorbox}[calloutstyle,title={Box~\thecallout\quad #2}]}
		{\end{tcolorbox}\end{figure}}

\newenvironment{calloutfloatplain}[1][tbp]% [placement]
{\refstepcounter{callout}%
	\begin{figure}[#1]\centering
		\begin{tcolorbox}[calloutstyle,title={Box~\thecallout}]}
		{\end{tcolorbox}\end{figure}}

% ----- Bibliography helper -----
% The class does NOT force a biblatex style to avoid coupling; users should load biblatex in their document and add resources.
% Example in your .tex:
%   \usepackage[backend=biber,style=alphabetic]{biblatex}
%   \addbibresource{references.bib}
%   ... \printbibliography

% ----- End of class -----
\endinput