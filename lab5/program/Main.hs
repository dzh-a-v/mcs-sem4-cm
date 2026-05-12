module Main where

import Text.Printf

-- ============================================================
-- Lab 6. Numerical integration
-- Variant 18
--
-- Function: y = arcctg(2^x)
-- Interval: [1; 3]
-- Newton-Cotes: n = 4
-- Gauss: n = 5
--
-- Console output uses Latin letters only.
-- ============================================================

a :: Double
a = 1.0

b :: Double
b = 3.0

-- For x in [1; 3], 2^x > 0.
-- arcctg(t) = atan(1 / t)
-- Therefore arcctg(2^x) = atan(1 / 2^x).
f :: Double -> Double
f x = atan (1.0 / (2.0 ** x))

-- ============================================================
-- Constants for error bounds.
--
-- M1 = max |f'(x)| on [1; 3]
-- M2 = max |f''(x)| on [1; 3]
-- M4 = max |f''''(x)| on [1; 3]
-- M6 = max |f''''''(x)| on [1; 3]
--
-- These constants are used only in theoretical error estimates
-- from the laboratory formulas.
-- ============================================================

m1 :: Double
m1 = 0.2772588722239781

m2 :: Double
m2 = 0.12011325347892068

m4 :: Double
m4 = 0.15733720319422967

m6 :: Double
m6 = 0.3128065157082597

-- ============================================================
-- Task 1. Simple composite quadrature formulas
-- ============================================================

leftRect :: Int -> (Double -> Double) -> Double -> Double -> Double
leftRect n func l r =
    let h = (r - l) / fromIntegral n
        xs = [l + fromIntegral i * h | i <- [0 .. n - 1]]
    in h * sum (map func xs)

rightRect :: Int -> (Double -> Double) -> Double -> Double -> Double
rightRect n func l r =
    let h = (r - l) / fromIntegral n
        xs = [l + fromIntegral i * h | i <- [1 .. n]]
    in h * sum (map func xs)

middleRect :: Int -> (Double -> Double) -> Double -> Double -> Double
middleRect n func l r =
    let h = (r - l) / fromIntegral n
        xs = [l + (fromIntegral i + 0.5) * h | i <- [0 .. n - 1]]
    in h * sum (map func xs)

trapezoid :: Int -> (Double -> Double) -> Double -> Double -> Double
trapezoid n func l r =
    let h = (r - l) / fromIntegral n
        innerSum = sum [func (l + fromIntegral i * h) | i <- [1 .. n - 1]]
    in h * ((func l + func r) / 2.0 + innerSum)

simpsonComposite :: Int -> (Double -> Double) -> Double -> Double -> Double
simpsonComposite nRaw func l r =
    let n = if even nRaw then nRaw else nRaw + 1
        h = (r - l) / fromIntegral n
        oddSum = sum [func (l + fromIntegral i * h) | i <- [1,3 .. n - 1]]
        evenSum = sum [func (l + fromIntegral i * h) | i <- [2,4 .. n - 2]]
    in h / 3.0 * (func l + func r + 4.0 * oddSum + 2.0 * evenSum)

leftRightErrorBound :: Int -> Double
leftRightErrorBound n =
    let h = (b - a) / fromIntegral n
    in (b - a) * h / 2.0 * m1

middleErrorBound :: Int -> Double
middleErrorBound n =
    let h = (b - a) / fromIntegral n
    in (b - a) * h * h / 24.0 * m2

trapezoidErrorBound :: Int -> Double
trapezoidErrorBound n =
    let h = (b - a) / fromIntegral n
    in (b - a) * h * h / 12.0 * m2

simpsonErrorBound :: Int -> Double
simpsonErrorBound nRaw =
    let n = if even nRaw then nRaw else nRaw + 1
        h = (b - a) / fromIntegral n
    in (b - a) * h ** 4.0 / 180.0 * m4

-- ============================================================
-- Task 2. Newton-Cotes formulas
--
-- General form:
-- Integral approx Bn * h * sum(ai * f(xi))
--
-- h = (b - a) / n
-- xi = a + i * h
-- ============================================================

newtonCotesData :: Int -> Maybe (Double, [Double], Double)
newtonCotesData n =
    case n of
        -- n = 1:
        -- B = 1/2, coefficients: 1, 1
        -- error bound: h^3 / 12 * M2
        1 -> Just (1.0 / 2.0, [1.0, 1.0], 1.0 / 12.0)

        -- n = 2:
        -- B = 1/3, coefficients: 1, 4, 1
        -- error bound: h^5 / 90 * M4
        2 -> Just (1.0 / 3.0, [1.0, 4.0, 1.0], 1.0 / 90.0)

        -- n = 3:
        -- B = 3/8, coefficients: 1, 3, 3, 1
        -- error bound: 3h^5 / 80 * M4
        3 -> Just (3.0 / 8.0, [1.0, 3.0, 3.0, 1.0], 3.0 / 80.0)

        -- n = 4:
        -- B = 2/45, coefficients: 7, 32, 12, 32, 7
        -- error bound: 8h^7 / 945 * M6
        4 -> Just (2.0 / 45.0, [7.0, 32.0, 12.0, 32.0, 7.0], 8.0 / 945.0)

        -- n = 5:
        -- B = 5/288, coefficients: 19, 75, 50, 50, 75, 19
        -- error bound: 275h^7 / 12096 * M6
        5 -> Just (5.0 / 288.0, [19.0, 75.0, 50.0, 50.0, 75.0, 19.0], 275.0 / 12096.0)

        _ -> Nothing

newtonCotes :: Int -> (Double -> Double) -> Double -> Double -> Double
newtonCotes n func l r =
    case newtonCotesData n of
        Nothing -> error "Unsupported Newton-Cotes n"
        Just (bn, coeffs, _) ->
            let h = (r - l) / fromIntegral n
                xs = [l + fromIntegral i * h | i <- [0 .. n]]
                values = map func xs
            in bn * h * sum (zipWith (*) coeffs values)

newtonCotesErrorBound :: Int -> Double
newtonCotesErrorBound n =
    case newtonCotesData n of
        Nothing -> error "Unsupported Newton-Cotes n"
        Just (_, _, c) ->
            let h = (b - a) / fromIntegral n
            in case n of
                1 -> c * h ** 3.0 * m2
                2 -> c * h ** 5.0 * m4
                3 -> c * h ** 5.0 * m4
                4 -> c * h ** 7.0 * m6
                5 -> c * h ** 7.0 * m6
                _ -> error "Unsupported Newton-Cotes n"

-- ============================================================
-- Task 3. Gauss formulas
--
-- Nodes and weights are taken from the Gauss table
-- for the interval [-1; 1].
--
-- For [a; b]:
-- ti = (a + b) / 2 + (b - a) / 2 * xi
--
-- Integral approx (b - a) / 2 * sum(Ai * f(ti))
-- ============================================================

gaussData :: Int -> Maybe [(Double, Double)]
gaussData n =
    case n of
        2 -> Just
            [ (-0.577350, 1.0)
            , ( 0.577350, 1.0)
            ]

        3 -> Just
            [ (-0.774597, 5.0 / 9.0)
            , ( 0.0,      8.0 / 9.0)
            , ( 0.774597, 5.0 / 9.0)
            ]

        4 -> Just
            [ (-0.861136, 0.347855)
            , (-0.339981, 0.652145)
            , ( 0.339981, 0.652145)
            , ( 0.861136, 0.347855)
            ]

        5 -> Just
            [ (-0.906180, 0.236927)
            , (-0.538469, 0.478629)
            , ( 0.0,      0.568889)
            , ( 0.538469, 0.478629)
            , ( 0.906180, 0.236927)
            ]

        _ -> Nothing

gauss :: Int -> (Double -> Double) -> Double -> Double -> Double
gauss n func l r =
    case gaussData n of
        Nothing -> error "Unsupported Gauss n"
        Just pairs ->
            let mid = (l + r) / 2.0
                half = (r - l) / 2.0
                oneTerm (x, weight) = weight * func (mid + half * x)
            in half * sum (map oneTerm pairs)

-- ============================================================
-- Printing
-- ============================================================

printLine :: IO ()
printLine = putStrLn "------------------------------------------------------------"

printSimpleResult :: String -> Int -> Double -> Double -> IO ()
printSimpleResult name n value bound = do
    printf "%-20s n = %-4d value = %.12f  error_bound = %.12e\n"
        name n value bound

printNCResult :: Int -> Double -> Double -> IO ()
printNCResult n value bound = do
    printf "Newton-Cotes         n = %-4d value = %.12f  error_bound = %.12e\n"
        n value bound

printGaussResult :: Int -> Double -> IO ()
printGaussResult n value = do
    printf "Gauss                n = %-4d value = %.12f\n"
        n value

-- ============================================================
-- Main
-- ============================================================

main :: IO ()
main = do
    putStrLn "Lab 5. Numerical integration"
    putStrLn "Variant 18"
    putStrLn "Function: y = arcctg(2^x)"
    putStrLn "Interval: [1; 3]"
    printLine

    putStrLn "Task 1. Simple composite formulas with different partition steps"
    putStrLn "Columns: formula, n, approximate value, theoretical error bound"
    printLine

    let parts = [4, 8, 16, 32]

    mapM_
        (\n -> do
            let value = leftRect n f a b
            let bound = leftRightErrorBound n
            printSimpleResult "Left rectangles" n value bound)
        parts

    printLine

    mapM_
        (\n -> do
            let value = rightRect n f a b
            let bound = leftRightErrorBound n
            printSimpleResult "Right rectangles" n value bound)
        parts

    printLine

    mapM_
        (\n -> do
            let value = middleRect n f a b
            let bound = middleErrorBound n
            printSimpleResult "Middle rectangles" n value bound)
        parts

    printLine

    mapM_
        (\n -> do
            let value = trapezoid n f a b
            let bound = trapezoidErrorBound n
            printSimpleResult "Trapezoid" n value bound)
        parts

    printLine

    mapM_
        (\n -> do
            let value = simpsonComposite n f a b
            let bound = simpsonErrorBound n
            printSimpleResult "Simpson" n value bound)
        parts

    printLine

    putStrLn "Task 2. Newton-Cotes formulas with different n"
    putStrLn "Columns: n, approximate value, theoretical error bound"
    printLine

    mapM_
        (\n -> do
            let value = newtonCotes n f a b
            let bound = newtonCotesErrorBound n
            printNCResult n value bound)
        [1 .. 5]

    printLine

    putStrLn "Task 3. Gauss formulas with different number of nodes"
    putStrLn "Columns: nodes count, approximate value"
    printLine

    mapM_
        (\n -> do
            let value = gauss n f a b
            printGaussResult n value)
        [2 .. 5]

    printLine

    putStrLn "Variant-specific results"
    let variantNC = newtonCotes 4 f a b
    let variantNCBound = newtonCotesErrorBound 4
    let variantGauss = gauss 5 f a b

    printf "Newton-Cotes n = 4 value = %.12f\n" variantNC
    printf "Newton-Cotes n = 4 error_bound = %.12e\n" variantNCBound
    printf "Gauss n = 5 value = %.12f\n" variantGauss

    printLine
    putStrLn "Done."