\name{housing}
\alias{housing}
\docType{data}
\title{
Housing data from kaggle
}
\description{
A subdata from kaggle "Get start" competition
}
\usage{data("housing")}
\format{
  A data frame with 10153 observations on the following 140 variables.
  \describe{
    \item{\code{floor}}{for apartments, floor of the building}
    \item{\code{area_m}}{Area, sq.m.}
    \item{\code{green_zone_part}}{Proportion of area of greenery in the total area}
    \item{\code{indust_part}}{Share of industrial zones in area of the total area}
    \item{\code{preschool_quota}}{Number of seats in pre-school organizations}
    \item{\code{preschool_education_centers_raion}}{Number of pre-school  institutions}
    \item{\code{school_quota}}{Number of high school seats in area}
    \item{\code{school_education_centers_raion}}{Number of  high school institutions}
    \item{\code{school_education_centers_top_20_raion}}{Number of high schools of the top 20 best schools in Moscow}
    \item{\code{healthcare_centers_raion}}{Number of healthcare centers in district	}
    \item{\code{university_top_20_raion}}{Number of higher education institutions in the top ten ranking of the Federal rank}
    \item{\code{sport_objects_raion}}{Number of higher education institutions}
    \item{\code{additional_education_raion}}{Number of additional education organizations}
    \item{\code{culture_objects_top_25_raion}}{Number of  objects of cultural heritage}
    \item{\code{shopping_centers_raion}}{Number of malls and shopping centres in district}
    \item{\code{office_raion}}{Number of malls and shopping centres in district}
    \item{\code{build_count_block}}{Share of block buildings}
    \item{\code{build_count_wood}}{Share of wood buildings}
    \item{\code{build_count_frame}}{Share of frame buildings}
    \item{\code{build_count_brick}}{Share of brick buildings}
    \item{\code{build_count_monolith}}{Share of monolith buildings}
    \item{\code{build_count_panel}}{Share of panel buildings}
    \item{\code{build_count_foam}}{Share of foam buildings}
    \item{\code{build_count_slag}}{Share of slag buildings}
    \item{\code{build_count_before_1920}}{Share of before_1920 buildings}
    \item{\code{build_count_1921.1945	}}{Share of 1921-1945 buildings}
    \item{\code{build_count_1946.1970}}{Share of 1946-1970 buildings}
    \item{\code{build_count_1971.1995}}{Share of 1971-1995 buildings}
    \item{\code{build_count_after_1995}}{Share of after_1995 buildings}
    \item{\code{kindergarten_km}}{Distance to kindergarten}
    \item{\code{school_km}}{Distance to high school}
    \item{\code{park_km}}{Distance to park}
    \item{\code{green_zone_km}}{Distance to green zone}
    \item{\code{industrial_km}}{Distance to industrial zone}
    \item{\code{water_treatment_km}}{Distance to water treatment}
    \item{\code{cemetery_km}}{Distance to the cemetery}
    \item{\code{incineration_km}}{Distance to the incineration}
    \item{\code{railroad_station_walk_min}}{Time to the railroad station (walk)}
    \item{\code{railroad_station_avto_km}}{Distance to the railroad station (avto)}
    \item{\code{railroad_station_avto_min}}{Time to the railroad station (avto)}
    \item{\code{public_transport_station_min_walk}}{Time to the public transport station (walk)}
    \item{\code{water_km}}{Distance to the water reservoir / river}
    \item{\code{mkad_km}}{Distance to MKAD (Moscow Circle Auto Road)}
    \item{\code{big_road1_km}}{Distance to Nearest major road}
    \item{\code{big_road2_km}}{The distance to next distant major road}
    \item{\code{railroad_km}}{Distance to the railway / Moscow Central Ring / open areas Underground}
    \item{\code{bus_terminal_avto_km}}{Distance to bus terminal (avto)}
    \item{\code{oil_chemistry_km}}{Distance to dirty industries}
    \item{\code{nuclear_reactor_km}}{Distance to nuclear reactor}
    \item{\code{radiation_km}}{Distance to burial of radioactive waste}
    \item{\code{power_transmission_line_km}}{Distance to power transmission line}
    \item{\code{thermal_power_plant_km}}{Distance to thermal power plant}
    \item{\code{ts_km}}{Distance to power station}
    \item{\code{big_market_km}}{Distance to grocery / wholesale markets}
    \item{\code{market_shop_km}}{Distance to markets and department stores}
    \item{\code{fitness_km}}{Distance to fitness}
    \item{\code{swim_pool_km}}{Distance to swimming pool}
    \item{\code{ice_rink_km}}{Distance to ice palace}
    \item{\code{stadium_km}}{Distance to stadium}
    \item{\code{basketball_km}}{Distance to the basketball courts}
    \item{\code{hospice_morgue_km}}{Distance to hospice/morgue}
    \item{\code{detention_facility_km}}{Distance to detention facility}
    \item{\code{public_healthcare_km}}{Distance to public healthcare}
    \item{\code{university_km}}{Distance to universities}
    \item{\code{workplaces_km}}{Distance to workplaces}
    \item{\code{shopping_centers_km}}{Distance to shopping centers}
    \item{\code{office_km}}{Distance to business centers/ offices}
    \item{\code{additional_education_km}}{Distance to additional education}
    \item{\code{preschool_km}}{Distance to preschool education organizations}
    \item{\code{big_church_km}}{Distance to large church}
    \item{\code{church_synagogue_km}}{Distance to Christian chirches and Synagogues}
    \item{\code{mosque_km}}{Distance to mosques}
    \item{\code{theater_km}}{Distance to theater}
    \item{\code{museum_km}}{Distance to museums}
    \item{\code{exhibition_km}}{Distance to exhibition}
    \item{\code{catering_km}}{Distance to catering}
    \item{\code{green_part_500}}{The share of green zones in 500 meters zone}
    \item{\code{prom_part_500}}{The share of industrial zones in 500 meters zone}
    \item{\code{office_count_500}}{The number of office space in 500 meters zone}
    \item{\code{office_sqm_500}}{The square of office space in 500 meters zone}
    \item{\code{trc_count_500}}{The number of shopping malls in 500 meters zone}
    \item{\code{trc_sqm_500}}{The square of shopping malls in 500 meters zone}
    \item{\code{cafe_count_500_na_price}}{Cafes and restaurant bill N/A in 500 meters zone}
    \item{\code{cafe_count_500_price_500}}{Cafes and restaurant bill, average under 500 in 500 meters zone}
    \item{\code{cafe_count_500_price_1000}}{Cafes and restaurant bill, average  500-1000 in 500 meters zone}
    \item{\code{cafe_count_500_price_1500}}{Cafes and restaurant bill, average  1000-1500 in 500 meters zone}
    \item{\code{cafe_count_500_price_2500}}{Cafes and restaurant bill, average  1500-2500 in 500 meters zone}
    \item{\code{cafe_count_500_price_4000}}{Cafes and restaurant bill, average  2500-4000 in 500 meters zone}
    \item{\code{cafe_count_500_price_high}}{Cafes and restaurant bill, average  over 4000 in 500 meters zone}
    \item{\code{big_church_count_500}}{The number of big churchs in 500 meters zone}
    \item{\code{church_count_500}}{The number of churchs in 500 meters zone}
    \item{\code{mosque_count_500}}{The number of mosques in 500 meters zone}
    \item{\code{leisure_count_500}}{The number of leisure facilities in 500 meters zone}
    \item{\code{sport_count_500}}{The number of sport facilities in 500 meters zone}
    \item{\code{market_count_500}}{The number of markets in 500 meters zone}
    \item{\code{green_part_1000}}{The share of green zones in 1000 meters zone}
    \item{\code{prom_part_1000}}{The share of industrial zones in 1000 meters zone}
    \item{\code{office_sqm_1000}}{The square of office space in 1000 meters zone}
    \item{\code{trc_count_1000}}{The number of shopping malls in 1000 meters zone}
    \item{\code{trc_sqm_1000}}{The square of shopping malls in 1000 meters zone}
    \item{\code{cafe_count_1000_na_price}}{Cafes and restaurant bill N/A in 1000 meters zone}
    \item{\code{cafe_count_1000_price_high}}{Cafes and restaurant bill, average  over 4000 in 1000 meters zone}
    \item{\code{big_church_count_1000}}{The number of big churchs in 1000 meters zone}
    \item{\code{mosque_count_1000}}{The number of mosques in 1000 meters zone}
    \item{\code{leisure_count_1000}}{The number of leisure facilities in 1000 meters zone}
    \item{\code{sport_count_1000}}{The number of sport facilities in 1000 meters zone}
    \item{\code{market_count_1000}}{The number of markets in 1000 meters zone}
    \item{\code{green_part_1500}}{The share of green zones in 1500 meters zone}
    \item{\code{prom_part_1500}}{The share of industrial zones in 1500 meters zone}
    \item{\code{office_sqm_1500}}{The square of office space in 1500 meters zone}
    \item{\code{trc_count_1500}}{	The number of shopping malls in 1500 meters zone}
    \item{\code{trc_sqm_1500}}{The square of shopping malls in 1500 meters zone}
    \item{\code{cafe_count_1500_price_high}}{Cafes and restaurant bill, average  over 4000 in 1500 meters zone	}
    \item{\code{mosque_count_1500}}{The number of mosques in 1500 meters zone}
    \item{\code{sport_count_1500}}{The number of sport facilities in 1500 meters zone}
    \item{\code{market_count_1500}}{The number of markets in 1500 meters zone}
    \item{\code{green_part_2000}}{The share of green zones in 2000 meters zone}
    \item{\code{prom_part_2000}}{The share of industrial zones in 2000 meters zone}
    \item{\code{office_sqm_2000}}{The square of office space in 2000 meters zone}
    \item{\code{trc_count_2000}}{The number of shopping malls in 2000 meters zone}
    \item{\code{trc_sqm_2000}}{The square of shopping malls in 2000 meters zone}
    \item{\code{mosque_count_2000}}{The number of mosques in 2000 meters zone}
    \item{\code{sport_count_2000}}{The number of sport facilities in 2000 meters zone}
    \item{\code{market_count_2000}}{The number of markets in 2000 meters zone}
    \item{\code{green_part_3000}}{The share of green zones in 3000 meters zone}
    \item{\code{prom_part_3000}}{The share of industrial zones in 3000 meters zone}
    \item{\code{office_sqm_3000}}{The square of office space in 3000 meters zone}
    \item{\code{trc_count_3000}}{The number of shopping malls in 3000 meters zone}
    \item{\code{trc_sqm_3000}}{The square of shopping malls in 3000 meters zone}
    \item{\code{mosque_count_3000}}{The number of mosques in 3000 meters zone}
    \item{\code{sport_count_3000}}{The number of sport facilities in 3000 meters zone}
    \item{\code{market_count_3000}}{The number of markets in 3000 meters zone}
    \item{\code{green_part_5000}}{The share of green zones in 5000 meters zone}
    \item{\code{prom_part_5000}}{The share of industrial zones in 5000 meters zone}
    \item{\code{trc_count_5000}}{The number of shopping malls in 5000 meters zone}
    \item{\code{trc_sqm_5000}}{The square of shopping malls in 5000 meters zone}
    \item{\code{mosque_count_5000}}{The number of mosques in 5000 meters zone}
    \item{\code{sport_count_5000}}{The number of sport facilities in 5000 meters zone}
    \item{\code{market_count_5000}}{The number of markets in 5000 meters zone}
    \item{\code{price_doc}}{I don't know}
    }
}
\source{
 www.kaggle.com
}
\examples{
data(housing)
}
\keyword{datasets}
