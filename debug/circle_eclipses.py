"""
There are n circles with known position and radius. Some circles overlap.
Calculate the visible relative area of each circle.
Complete the function visible_percentages().
Consider the funtion intersection() to be correctly implemented.
"""

class TCircle():

    def __init__(self, x, y, r, height):
        self.x = x  # x-coordinate
        self.y = y  # y-coordinate
        self.r = r  # radius
        self.height = height  #  determines visibility precedence (higher circles obscure lower ones)


def intersection(eclipsing_circle, eclipsed_circle):
    """ Returns the relative area of eclipsed_circle which is eclipsed by eclipsing_circle."""
    # Placeholder. Consider this funtion to be correctly implemented
    # For demonstration, let's make it return a simple value based on some overlap idea.
    # In a real scenario, this would involve geometric calculations.
    # A simple placeholder: if they are close and eclipsing_circle is larger,
    # it might eclipse more. This is NOT a real geometric calculation.

    # Example of a *very simplified* placeholder behavior for illustration:
    # If centers are close and eclipsing_circle is larger, assume some overlap.
    dist_sq = (eclipsing_circle.x - eclipsed_circle.x)**2 + \
              (eclipsing_circle.y - eclipsed_circle.y)**2

    # If the circles are very far apart, no intersection
    if dist_sq > (eclipsing_circle.r + eclipsed_circle.r)**2:
        return 0.0

    # If one circle is completely inside the other and eclipsing is larger
    if dist_sq < (eclipsing_circle.r - eclipsed_circle.r)**2:
        if eclipsing_circle.r >= eclipsed_circle.r:
            return 1.0 # eclipsed_circle is fully covered
        else:
            return 0.0 # eclipsed_circle covers eclipsing_circle, but not vice-versa

    # Otherwise, assume some partial overlap.
    # This is where a real geometric intersection calculation would be.
    # For this example, let's return a dummy value based on radius ratio as a placeholder.
    # This is NOT mathematically sound for actual intersection area.
    if eclipsing_circle.r > eclipsed_circle.r:
        return min(1.0, 0.5 * (eclipsing_circle.r / eclipsed_circle.r)) # Dummy value
    else:
        return min(1.0, 0.2 * (eclipsing_circle.r / eclipsed_circle.r)) # Dummy value



def visible_percentages_gemini(circle_list):
    """Returns a list of percentages for the visible relative area of each circle."""

    # Create a list of (initial_index, TCircle) tuples to preserve original order
    indexed_circles = [(i, circle) for i, circle in enumerate(circle_list)]

    # Sort circles by height in ascending order
    # Circles with higher height will be considered later and can obscure earlier ones
    indexed_circles.sort(key=lambda x: x[1].height)

    # Initialize visible areas for all circles to 1.0 (100%)
    # This list will store results in the original order
    visible_areas = [1.0] * len(circle_list)

    # Iterate through the sorted circles
    for i, (original_idx_i, circle_i) in enumerate(indexed_circles):
        # If a circle is already completely obscured by prior (higher) circles, skip it.
        # This logic needs careful consideration. If we subtract directly, we don't need this check here.
        # We'll subtract and ensure we don't go below zero.

        current_circle_visible_area = 1.0

        # Check for eclipsing by circles with higher height
        for j, (original_idx_j, circle_j) in enumerate(indexed_circles):
            if circle_j.height > circle_i.height:
                # Calculate the proportion of circle_i eclipsed by circle_j
                # IMPORTANT: This assumes 'intersection' returns the *additional*
                # eclipsed area by circle_j, not already covered by other higher circles.
                # If 'intersection' returns total overlap area, this will overcount.
                eclipsed_proportion = intersection(circle_j, circle_i)
                current_circle_visible_area -= eclipsed_proportion

                # Ensure the visible area doesn't drop below zero
                current_circle_visible_area = max(0.0, current_circle_visible_area)

        # Store the final visible area for circle_i in its original index
        visible_areas[original_idx_i] = current_circle_visible_area

    return visible_areas

def visible_percentages_perplexity(circle_list):
    """Returns a list of percentages for the visible relative area of each circle."""

    # Create a list of (initial_index, TCircle) tuples to preserve original order
    indexed_circles = [(i, circle) for i, circle in enumerate(circle_list)]

    # Sort circles by height in ascending order
    indexed_circles.sort(key=lambda x: x[1].height)

    # Initialize visible areas for all circles to 1.0 (100%)
    visible_areas = [1.0] * len(circle_list)

    # Iterate through the sorted circles
    for i, (original_idx_i, circle_i) in enumerate(indexed_circles):
        current_circle_visible_area = 1.0

        # Track the cumulative eclipsed area
        cumulative_eclipsed_area = 0.0

        # Check for eclipsing by circles with higher height
        for j, (original_idx_j, circle_j) in enumerate(indexed_circles):
            if circle_j.height > circle_i.height:
                # Calculate the proportion of circle_i eclipsed by circle_j
                eclipsed_proportion = intersection(circle_j, circle_i)

                # Subtract only the additional eclipsed area not already covered
                additional_eclipsed_area = max(0.0, eclipsed_proportion - cumulative_eclipsed_area)
                cumulative_eclipsed_area += additional_eclipsed_area

                # Ensure the visible area doesn't drop below zero
                current_circle_visible_area -= additional_eclipsed_area
                current_circle_visible_area = max(0.0, current_circle_visible_area)

        # Store the final visible area for circle_i in its original index
        visible_areas[original_idx_i] = current_circle_visible_area

    return visible_areas

# Python
def visible_percentages_copilot(circle_list):
    """Returns a list of percentages for the visible relative area of each circle."""

    # Create a list of (initial_index, TCircle) tuples to preserve original order
    indexed_circles = [(i, circle) for i, circle in enumerate(circle_list)]

    # Sort circles by height in ascending order
    indexed_circles.sort(key=lambda x: x[1].height)

    # Initialize visible areas for all circles to 1.0 (100%)
    visible_areas = [1.0] * len(circle_list)

    # Iterate through the sorted circles
    for i, (original_idx_i, circle_i) in enumerate(indexed_circles):
        current_circle_visible_area = 1.0

        # Track the cumulative eclipsed area
        cumulative_eclipsed_area = 0.0

        # Check for eclipsing by circles with higher height
        for j, (original_idx_j, circle_j) in enumerate(indexed_circles):
            if circle_j.height > circle_i.height:
                # Calculate the proportion of circle_i eclipsed by circle_j
                eclipsed_proportion = intersection(circle_j, circle_i)

                # Subtract only the additional eclipsed area not already covered
                additional_eclipsed_area = max(0.0, eclipsed_proportion - cumulative_eclipsed_area)
                cumulative_eclipsed_area += additional_eclipsed_area

                # Ensure the visible area doesn't drop below zero
                current_circle_visible_area -= additional_eclipsed_area
                current_circle_visible_area = max(0.0, current_circle_visible_area)

        # Store the final visible area for circle_i in its original index
        visible_areas[original_idx_i] = current_circle_visible_area

    return visible_areas


# Example circles
circle1 = TCircle(10, 20, 25, 1)
circle2 = TCircle(-10, 0, 35, 2)
circle3 = TCircle(10, 25, 15, 3)
circle4 = TCircle(12, 28, 22, 4)
circle5 = TCircle(12.1, 28.1, 2.2, 5)

circle_list = [circle3, circle5, circle4, circle2, circle1]  # should be something like [.42, .68, .55, 1.0]

print(visible_percentages_gemini(circle_list))
print(visible_percentages_perplexity(circle_list))
print(visible_percentages_copilot(circle_list))
