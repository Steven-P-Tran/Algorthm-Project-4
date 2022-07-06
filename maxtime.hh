///////////////////////////////////////////////////////////////////////////////
// maxtime.hh
//
// Compute the set of rides that maximizes the time spent at rides, within a given budget
// with the dynamic method or exhaustive search.
//
///////////////////////////////////////////////////////////////////////////////


#pragma once


#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <queue>
#include <sstream>
#include <string>
#include <vector>


// One ride item available for purchase.
class RideItem
{
	//
	public:

		//
		RideItem
		(
			const std::string& description,
			size_t cost_dollars,
			double time_minutes
		)
			:
			_description(description),
			_cost_dollars(cost_dollars),
			_time_minutes(time_minutes)
		{
			assert(!description.empty());
			assert(cost_dollars > 0);
		}

		//
		const std::string& description() const { return _description; }
		int cost() const { return _cost_dollars; }
		double defense() const { return _time_minutes; }

	//
	private:

		// Human-readable description of the ride, e.g. "new enchanted world". Must be non-empty.
		std::string _description;

		// Ride cost, in units of dollars; Must be positive
		int _cost_dollars;

		// Ride time in minutes; most be non-negative.
		double _time_minutes;
};


// Alias for a vector of shared pointers to RideItem objects.
typedef std::vector<std::shared_ptr<RideItem>> RideVector;


// Load all the valid ride items from the CSV database
// ride items that are missing fields, or have invalid values, are skipped.
// Returns nullptr on I/O error.
std::unique_ptr<RideVector> load_ride_database(const std::string& path)
{
	std::unique_ptr<RideVector> failure(nullptr);

	std::ifstream f(path);
	if (!f)
	{
		std::cout << "Failed to load ride database; Cannot open file: " << path << std::endl;
		return failure;
	}

	std::unique_ptr<RideVector> result(new RideVector);

	size_t line_number = 0;
	for (std::string line; std::getline(f, line); )
	{
		line_number++;

		// First line is a header row
		if ( line_number == 1 )
		{
			continue;
		}

		std::vector<std::string> fields;
		std::stringstream ss(line);

		for (std::string field; std::getline(ss, field, '^'); )
		{
			fields.push_back(field);
		}

		if (fields.size() != 3)
		{
			std::cout
				<< "Failed to load ride database: Invalid field count at line " << line_number << "; Want 3 but got " << fields.size() << std::endl
				<< "Line: " << line << std::endl
				;
			return failure;
		}

		std::string
			descr_field = fields[0],
			cost_dollars_field = fields[1],
			time_minutes_field = fields[2]
			;

		auto parse_dbl = [](const std::string& field, double& output)
		{
			std::stringstream ss(field);
			if ( ! ss )
			{
				return false;
			}

			ss >> output;

			return true;
		};

		std::string description(descr_field);
		double cost_dollars, time_minutes;
		if (
			parse_dbl(cost_dollars_field, cost_dollars)
			&& parse_dbl(time_minutes_field, time_minutes)
		)
		{
			result->push_back(
				std::shared_ptr<RideItem>(
					new RideItem(
						description,
						cost_dollars,
						time_minutes
					)
				)
			);
		}
	}

	f.close();

	return result;
}


// Convenience function to compute the total cost and time in a RideVector.
// Provide the RideVector as the first argument
// The next two arguments will return the cost and time back to the caller.
void sum_ride_vector
(
	const RideVector& rides,
	int& total_cost,
	double& total_time
)
{
	total_cost = total_time = 0;
	for (auto& ride : rides)
	{
		total_cost += ride->cost();
		total_time += ride->defense();
	}
}


// Convenience function to print out each RideItem in a RideVector,
// followed by the total kilocalories and protein in it.
void print_ride_vector(const RideVector& rides)
{
	std::cout << "*** ride Vector ***" << std::endl;

	if ( rides.size() == 0 )
	{
		std::cout << "[empty ride list]" << std::endl;
	}
	else
	{
		for (auto& ride : rides)
		{
			std::cout
				<< "Ye olde " << ride->description()
				<< " ==> "
				<< "Cost of " << ride->cost() << " dollars"
				<< "; time points = " << ride->defense()
				<< std::endl
				;
		}

		int total_cost;
		double total_time;
		sum_ride_vector(rides, total_cost, total_time);
		std::cout
			<< "> Grand total cost: " << total_cost << " dollars" << std::endl
			<< "> Grand total defense: " << total_time
			<< std::endl
			;
	}
}


// Convenience function to print out a 2D cache, composed of a std::vector<std::vector<double>>
// For sanity, will refuse to print a cache that is too large.
// Hint: When running this program, you can redirect stdout to a file,
//	which may be easier to view and inspect than a terminal
void print_2d_cache(const std::vector<std::vector<double>>& cache)
{
	std::cout << "*** 2D Cache ***" << std::endl;

	if ( cache.size() == 0 )
	{
		std::cout << "[empty]" << std::endl;
	}
	else if ( cache.size() > 250 || cache[1].size() > 250 )
	{
		std::cout << "[too large]" << std::endl;
	}
	else
	{
		for ( const std::vector<double>& row : cache)
		{
			for ( double value : row )
			{
				std::cout << std::setw(5) << value;
			}
			std::cout << std::endl;
		}
	}
}

// Filter the vector source, i.e. create and return a new RideVector containing the subset of 
// the ride items in source that match given criteria.
// This is intended to:
//	1) filter out ride with zero or negative time that are irrelevant to our optimization
//	2) limit the size of inputs to the exhaustive search algorithm since it will probably be slow.
//
// Each ride item that is included must have at minimum min_time and at most max_time.
//	(i.e., each included ride item's time must be between min_time and max_time (inclusive).
//
// In addition, the the vector includes only the first total_size ride items that match these criteria.
std::unique_ptr<RideVector> filter_ride_vector
(
	const RideVector &source,
	double min_time,
	double max_time,
	int total_size
)
{
	std::unique_ptr<RideVector> subset(new RideVector());

	int total_rides_added = 0;
	for (unsigned long i = 0; total_rides_added < total_size && i < source.size(); i++) {
		std::shared_ptr<RideItem> current_ride = source[i];
		double current_ride_time = current_ride->defense();

		if (!(current_ride_time > 0 && current_ride_time >= min_time && current_ride_time <= max_time)) {
			continue;
		}
		subset->push_back(current_ride);
		total_rides_added++;
	}
	return subset;
}

void recursive_solver(
	const RideVector &rides,
	std::vector<std::vector<double>> &cache,
	RideVector &resultSet,
	int cost,
	int n
) {
	// If no more items are to be checked, then return
	if (n == 0) {
		return;
	}

	// If the current value and the value above are the same then the item wasn't added and we move back to the next item
	if (cache[n][cost] == cache[n - 1][cost]) {
		recursive_solver(rides, cache, resultSet, cost, n - 1);
	}
	// If the current item was added, then we add it to the resultSet, and move up the cache (lower idx) and remove the cost of the currently added Item from the cost available to be used.
	else {
		std::shared_ptr<RideItem> itemToAdd = rides[n - 1];
		resultSet.push_back(itemToAdd);

		recursive_solver(rides, cache, resultSet, cost - itemToAdd->cost(), n - 1);
	}

}

// Compute the optimal set of ride items with a dynamic algorithm.
// Specifically, among the ride items that fit within a total_cost budget,
// choose the selection of rides whose time is greatest.
// Repeat until no more ride items can be chosen, either because we've run out of ride items,
// or run out of dollars.
std::unique_ptr<RideVector> dynamic_max_time
(
	const RideVector& rides,
	int total_cost
)
{
	size_t total_rides = rides.size();

	// Creates and initializing the cache with zeros
	std::vector<std::vector<double>> cache;
	for (int i = 0, n = total_rides + 1; i < n; i++) {
		std::vector<double> subVector;
		for (int j = 0; j < total_cost + 1; j++) {
			subVector.push_back(0);
		}
		cache.push_back(subVector);
	}


	// Filling the cache
	for (int i = 1; i < cache.size(); i++) {
		std::shared_ptr<RideItem> currItem = rides[i - 1];
		for (int j = 0; j < total_cost + 1; j++) {
			double value = 0;

			if (j >= currItem->cost()) {
				value = cache[i - 1][j - currItem->cost()] + currItem->defense();
			}
			cache[i][j] = std::max(value, cache[i - 1][j]);
		}
	}

	// Solving the for best RideItems with the help of cache
	RideVector resultSet;
	recursive_solver(rides, cache, resultSet, total_cost, total_rides);

	return std::unique_ptr<RideVector>(new RideVector(resultSet));
}
std::vector<std::vector<RideItem>> getDefenseSubsets(std::vector<RideItem> source)
{
    std::vector<std::vector<RideItem>> subset, subTemp;

    std::vector<RideItem> temp;
    subset.push_back(temp);

    for(int i = 0; i < source.size(); i++)
    {
        subTemp = subset;
        for(int j = 0; j < subTemp.size(); j++){subTemp[j].push_back(source[i]);}
        for(int k = 0; k < subTemp.size(); k++){subset.push_back(subTemp[k]);}
    }
    return subset;
}

// Compute the optimal set of ride items with a exhaustive search algorithm.
// Specifically, among all subsets of ride items,
// return the subset whose dollars cost fits within the total_cost budget,
// and whose total time is greatest.
// To avoid overflow, the size of the ride items vector must be less than 64.
std::unique_ptr<RideVector> exhaustive_max_time
(
	const RideVector& rides,
	double total_cost
)
{
	size_t n = rides.size();
	// Returns nullptr if the rides' count >= 64
	if (n >= 64) { return nullptr; }

	std::unique_ptr<RideVector> best(nullptr);

	for (uint64_t bits = 0; bits <= (pow(2, n) - 1); bits++) {
		std::unique_ptr<RideVector> candidate(new RideVector());
		for (uint64_t j = 0; j <= (n - 1); j++) {
			if (((bits >> j) & 1) == 1) {
				candidate->push_back(rides[j]);
			}
		}

		int candidate_tc, best_tc;
		double candidate_tt, best_tt;
		sum_ride_vector(*candidate, candidate_tc, candidate_tt);

		if (best != nullptr) {
			sum_ride_vector(*best, best_tc, best_tt);
		}

		if (candidate_tc <= total_cost) {
			if (best == nullptr || (candidate_tt > best_tt)) {
				best = std::move(candidate);
			}
		}
	}

	return best;
}








