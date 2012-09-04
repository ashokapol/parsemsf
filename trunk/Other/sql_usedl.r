SELECT paam.PeptideID, p.Sequence, p.ProcessingNodeNumber AS SearchNode,
group_concat(paam.Position) AS Position,
group_concat(aam.ModificationName) AS ModificationName,
p.ConfidenceLevel
FROM PeptidesAminoAcidModifications paam, AminoAcidModifications aam, Peptides p
WHERE paam.AminoAcidModificationID = aam.AminoAcidModificationID AND
paam.PeptideID = p.PeptideID AND
p.ConfidenceLevel >=3
GROUP BY paam.PeptideID
ORDER BY paam.PeptideID;



SELECT p.PeptideID, p.Sequence, sum(ev.Area) AS Area
FROM Peptides p, Events ev, EventAreaAnnotations eva, PrecursorIonAreaSearchSpectra pis
WHERE ev.EventID = eva.EventID
AND eva.QuanResultID = pis.QuanResultID
AND pis.SearchSpectrumID = p.SpectrumID AND
p.ConfidenceLevel >=3
GROUP BY p.PeptideID;


SELECT p.PeptideID, p.Sequence, piqr.QuanChannelID, piqr.Area AS Area
FROM Peptides p, PrecursorIonAreaSearchSpectra pis,
PrecursorIonQuanResults piqr
WHERE piqr.QuanResultID = pis.QuanResultID
AND pis.SearchSpectrumID = p.SpectrumID AND
p.ConfidenceLevel >=1;


SELECT p.PeptideID, p.Sequence, piqr.QuanChannelID, p.ConfidenceLevel, piqr.Area AS Area
FROM Peptides p, PrecursorIonAreaSearchSpectra pis,
PrecursorIonQuanResults piqr
WHERE piqr.QuanResultID = pis.QuanResultID
AND pis.SearchSpectrumID = p.SpectrumID;

SELECT p.PeptideID, p.Sequence, sum(ev.Area) AS Area
FROM Peptides p, Events ev, EventAreaAnnotations eva, PrecursorIonAreaSearchSpectra pis
WHERE ev.EventID = eva.EventID
AND eva.QuanResultID = pis.QuanResultID
AND pis.SearchSpectrumID = p.SpectrumID AND
p.ConfidenceLevel >=1;

SELECT p.PeptideID, p.Sequence, pias.QuanResultID, e.Area, piqr.QuanChannelID
FROM Peptides p, PrecursorIonAreaSearchSpectra pias, EventAreaAnnotations eaa, Events e,
PrecursorIonQuanResults piqr
WHERE pias.SearchSpectrumID = p.SpectrumID AND
eaa.QuanResultID = pias.QuanResultID AND
e.EventID = eaa.EventID AND
piqr.QuanResultID = pias.QuanResultID AND
p.ConfidenceLevel >=1;

SELECT p.PeptideID, p.Sequence, p.ConfidenceLevel, eva.IsotopePatternID, eva.QuanChannelID, e.Area
FROM Peptides p, PrecursorIonQuanResultsSearchSpectra piqrss,
EventAnnotations eva, Events e
WHERE  piqrss.SearchSpectrumID = p.SpectrumID AND
eva.QuanResultID = piqrss.QuanResultID AND
e.EventID = eva.EventID AND
p.ConfidenceLevel >=2;
--GROUP BY eva.IsotopePatternID;

