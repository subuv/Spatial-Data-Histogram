CREATE TABLE histogram_input(num_particles, bucket_width) AS VALUES (10000, 500);

CREATE OR REPLACE FUNCTION SDH(INT, INT) RETURNS INT AS 'spatialDataHistogram.so', 'spatialDataHistogram' LANGUAGE C STRICT IMMUTABLE;

SELECT SDH(num_particles, bucket_width) FROM histogram_input;