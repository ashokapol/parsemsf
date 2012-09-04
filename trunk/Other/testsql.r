select Peptides.PeptideID, ReporterIonQuanResults.Mass, ReporterIonQuanResults.Height
from Peptides, ReporterIonQuanResults, ReporterIonQuanResultsSearchSpectra
where Peptides.SpectrumID = ReporterIonQuanResultsSearchSpectra.SearchSpectrumID AND
ReporterIonQuanResultsSearchSpectra.SpectrumID = ReporterIonQuanResults.SpectrumID;

select Peptides.PeptideID, Peptides.Sequence, SpectrumHeaders.ScanEventID, ReporterIonQuanResults.Mass, ReporterIonQuanResults.Height
from Peptides, ReporterIonQuanResults, ReporterIonQuanResultsSearchSpectra, SpectrumHeaders
where Peptides.SpectrumID = ReporterIonQuanResultsSearchSpectra.SearchSpectrumID AND Peptides.SpectrumID=SpectrumHeaders.SpectrumID AND
ReporterIonQuanResultsSearchSpectra.SpectrumID = ReporterIonQuanResults.SpectrumID AND
Peptides.Sequence='GPSCQDCDTGYTR' AND Peptides.ConfidenceLevel=3;