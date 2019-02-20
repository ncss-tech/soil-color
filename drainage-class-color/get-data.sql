

SET search_path to osd, soilweb, public;

CREATE TEMP TABLE mlra_drainage_colors AS
SELECT mlra, membership, series, drainagecl, hzname, top, bottom, matrix_dry_color_hue, matrix_dry_color_value, matrix_dry_color_chroma, matrix_wet_color_hue, matrix_wet_color_value, matrix_wet_color_chroma
FROM osd_site JOIN osd_colors USING (series)
JOIN (
SELECT DISTINCT ON (series) series, mlra, membership
FROM mlra_overlap
ORDER BY series, membership DESC
) as a USING (series) 
ORDER BY mlra, series, top ASC ;

\copy mlra_drainage_colors TO 'mlra_drainage_colors.csv' CSV HEADER
