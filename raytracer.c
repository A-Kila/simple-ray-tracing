#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define WIDTH 1920
#define HEIGHT 1080

#define EPSILON 0.001f

#define DEBUG_IF(condition, message, ...) \
    if (condition) { \
        printf("DEBUG: %s--%d: ", __FILE__, __LINE__); \
        printf(message, ##__VA_ARGS__); \
        printf("\n"); \
    }

typedef struct vec3f_t vec3f_t;
typedef struct array_t array_t;
vec3f_t get_pixel_color(
    const vec3f_t *ray_origin, 
    const vec3f_t *ray_direction, 
    const array_t *spheres,
    const array_t *lights,
    const int depth);

// ------------------------------------------------------------------

struct vec3f_t
{
    float x; // r
    float y; // g
    float z; // b
};

vec3f_t vec3f_add(const vec3f_t a, const vec3f_t b)
{
    return (vec3f_t){ a.x + b.x, a.y + b.y, a.z + b.z };
}

vec3f_t vec3f_sub(const vec3f_t a, const vec3f_t b)
{
    return (vec3f_t){ a.x - b.x, a.y - b.y, a.z - b.z };
}

vec3f_t vec3f_scale(const vec3f_t v, float scalar)
{
    return (vec3f_t){ v.x * scalar, v.y * scalar, v.z * scalar };
}

vec3f_t vec3f_normalize(const vec3f_t v)
{
    float length = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
    
    if (length == 0) return (vec3f_t){ 0, 0, 0 };

    return (vec3f_t){ v.x / length, v.y / length, v.z / length };
}

float vec3f_dot(const vec3f_t a, const vec3f_t b)
{
    return a.x * b.x+ a.y * b.y + a.z * b.z;
}

float vec3f_length(const vec3f_t v)
{
    return sqrtf(v.x * v.x+ v.y * v.y + v.z * v.z);
}

// ------------------------------------------------------------------

struct array_t
{
    void *data;
    size_t size;
};

// ------------------------------------------------------------------

typedef struct
{
    vec3f_t color;
    float shininess;
    float refraction_index;
    float diffuse_constant;
    float specular_constant;
    float reflection_constant;
    float refraction_constant;
} material_t;

typedef struct
{
    vec3f_t position;
    float intensity;    
} light_t;

typedef struct
{
    vec3f_t center;
    float radius;
    material_t material;
} sphere_t;

bool sphere_is_ray_intersecting(
    const sphere_t *sphere, 
    const vec3f_t *ray_origin, 
    const vec3f_t *ray_direction, 
    vec3f_t *intersection_point)
{
    const vec3f_t ray_to_center_vec = vec3f_sub(sphere->center, *ray_origin);
    const float rc_dot = vec3f_dot(ray_to_center_vec, *ray_direction);

    if (rc_dot < 0) return false; // ray points away from the sphere
    if (vec3f_length(ray_to_center_vec) < sphere->radius) return false; // no need to render when we are inside the sphere

    const vec3f_t projection = vec3f_scale(*ray_direction, rc_dot);
    const vec3f_t projection_end = vec3f_add(*ray_origin, projection);
    const vec3f_t sphere_to_projection = vec3f_sub(sphere->center, projection_end);
    const float distance_squared = vec3f_dot(sphere_to_projection, sphere_to_projection);

    if (distance_squared > sphere->radius * sphere->radius) return false; // no intersection

    const float dist_projection_to_intersection = sqrtf(sphere->radius * sphere->radius - distance_squared);
    const float intersection_distance = rc_dot - dist_projection_to_intersection;
    const vec3f_t intersection_vector = vec3f_scale(*ray_direction, intersection_distance);
    *intersection_point = vec3f_add(*ray_origin, intersection_vector);

    return true;   
}

// ------------------------------------------------------------------

vec3f_t get_light_direction(const light_t *light, const vec3f_t *hit_point)
{
    vec3f_t hit_to_light = vec3f_sub(light->position, *hit_point);
    return vec3f_normalize(hit_to_light);

}

vec3f_t reflect_ray(const vec3f_t *incident, const vec3f_t *normal)
{
    float dot_product = vec3f_dot(*incident, *normal);
    vec3f_t scaled_normal = vec3f_scale(*normal, 2 * dot_product);
    return vec3f_sub(*incident, scaled_normal);
}

vec3f_t clamp_color(const vec3f_t *color)
{
    float max_color_value = fmaxf(fmaxf(color->x, color->y), color->z);

    if (max_color_value <= 0.f) return (vec3f_t){0.f, 0.f, 0.f};
    
    if (max_color_value > 1.f)
        return (vec3f_t){
            .x = color->x / max_color_value,
            .y = color->y / max_color_value,
            .z = color->z / max_color_value
        };

    return (vec3f_t){color->x, color->y, color->z};
}

bool cast_ray_to_spheres(const vec3f_t *ray_origin, 
                          const vec3f_t *ray_direction, 
                          const array_t *spheres,
                          vec3f_t *hit_point,
                          vec3f_t *normal,
                          material_t *out_material)
{
    static const float MAX_DISTANCE = 500.f;

    int distance = MAX_DISTANCE + 1;
    bool hit = false;

    for (int i = 0; i < spheres->size; i++)
    {
        const sphere_t *sphere = (sphere_t *)spheres->data + i;
        
        const vec3f_t ray_to_intersection = vec3f_sub(sphere->center, *ray_origin);
        if (vec3f_length(ray_to_intersection) >= distance) continue;

        vec3f_t intersection_point;
        if (sphere_is_ray_intersecting(sphere, ray_origin, ray_direction, &intersection_point))
        {
            distance = vec3f_length(intersection_point);
            *hit_point = intersection_point;
            
            *normal = vec3f_normalize(vec3f_sub(intersection_point, sphere->center));

            *out_material = sphere->material;
            hit = true;
        }
    }

    return hit;
}

bool is_shadow(
    const light_t *light,
    const vec3f_t *hit_point,
    const vec3f_t *light_dir,
    const vec3f_t *normal,
    const array_t *spheres)
{
    const float light_distance = vec3f_length(vec3f_sub(light->position, *hit_point));

    const float dot_product = vec3f_dot(*light_dir, *normal);
    const vec3f_t small_step = vec3f_scale(*normal, EPSILON); // to avoid self-intersection
    vec3f_t shadow_origin = dot_product < 0 ? vec3f_sub(*hit_point, small_step) : vec3f_add(*hit_point, small_step);

    vec3f_t shadow_hit_point, shadow_normal;
    material_t shadow_material;

    bool result = cast_ray_to_spheres(&shadow_origin, light_dir, spheres, &shadow_hit_point, &shadow_normal, &shadow_material);
    
    result = result && (vec3f_length(vec3f_sub(shadow_hit_point, shadow_origin)) < light_distance);
    
    return result;
}

float calculate_diffuse_lighting(
    const array_t *lights, 
    const vec3f_t *hit_point, 
    const vec3f_t *normal, 
    const material_t *material, 
    const array_t *spheres)
{   
    float diffuse_intensity = 0.f;

    for (int i = 0; i < lights->size; i++)
    {
        const vec3f_t light_dir = get_light_direction((light_t *)lights->data + i, hit_point);
        if (is_shadow((light_t *)lights->data + i, hit_point, &light_dir, normal, spheres)) continue;

        diffuse_intensity += ((light_t *)lights->data)[i].intensity * fmaxf(0.f, vec3f_dot(light_dir, *normal));
    }

    return (diffuse_intensity * material->diffuse_constant);
}

vec3f_t calculate_specular_lighting(
    const array_t *lights,
    const vec3f_t *ray_direction,
    const vec3f_t *hit_point,
    const vec3f_t *normal,
    const material_t *material,
    const array_t *spheres)
{
    const vec3f_t specular_color = {1.f, 1.f, 1.f};
    float specular_intensity = 0.f;

    for (int i = 0; i < lights->size; i++)
    {
        const light_t *light = (light_t *)lights->data + i;
        const vec3f_t light_dir = get_light_direction(light, hit_point);
        
        if (is_shadow((light_t *)lights->data + i, hit_point, &light_dir, normal, spheres)) continue;

        const vec3f_t reflected_ray = reflect_ray(&light_dir, normal);
        float specular_factor = fmaxf(0.f, vec3f_dot(*ray_direction, reflected_ray));

        specular_intensity += light->intensity * powf(specular_factor, material->shininess);
    }

    return vec3f_scale(specular_color, specular_intensity * material->specular_constant);
}

vec3f_t calculate_reflected_color(
    const vec3f_t *ray_direction,
    const vec3f_t *normal,
    const vec3f_t *hit_point,
    const array_t *spheres,
    const array_t *lights,
    const material_t *material,
    const int depth
)
{
    const vec3f_t reflected_ray = vec3f_normalize(reflect_ray(ray_direction, normal));
    const vec3f_t small_step = vec3f_scale(*normal, EPSILON); // to avoid self-intersection
    const float dot = vec3f_dot(reflected_ray, *normal);
    const vec3f_t ray_origin = dot < 0 ? vec3f_sub(*hit_point, small_step) : vec3f_add(*hit_point, small_step);

    const vec3f_t ref_color = get_pixel_color(&ray_origin, &reflected_ray, spheres, lights, depth + 1);

    return vec3f_scale(ref_color, material->reflection_constant);
}

vec3f_t refract_ray(const vec3f_t *ray_direction, const vec3f_t *normal, const float refraction_index)
{
    float cos_t1 = vec3f_dot(*ray_direction, *normal);
    float n1 = 1.f; // Air refraction index
    float n2 = refraction_index;

    vec3f_t norm = *normal;
    // if ray is coming out of the material, invert values
    if (cos_t1 < 0)
    {
        cos_t1 = -cos_t1;
        n1 = refraction_index;
        n2 = 1.f;
        norm = vec3f_scale(*normal, -1.f);
    }

    const float ratio = n1 / n2;
    const float cos_t2_squared = 1.f - ratio * ratio * (1.f - cos_t1 * cos_t1);

    if (cos_t2_squared < 0.f) return (vec3f_t){0.f, 0.f, 0.f};

    const float cos_t2 = sqrtf(cos_t2_squared);

    return vec3f_add(vec3f_scale(*ray_direction, ratio), vec3f_scale(norm, ratio * cos_t1 - cos_t2));   
}

vec3f_t calculate_refracted_color
(
    const vec3f_t *ray_direction,
    const vec3f_t *normal,
    const vec3f_t *hit_point,
    const array_t *spheres,
    const array_t *lights,
    const material_t *material,
    const int depth
)
{
    const vec3f_t refracted_ray = vec3f_normalize(refract_ray(ray_direction, normal, material->refraction_index));
    const vec3f_t small_step = vec3f_scale(*normal, EPSILON); // to avoid self-intersection
    const float dot = vec3f_dot(refracted_ray, *normal);
    const vec3f_t ray_origin = dot < 0 ? vec3f_sub(*hit_point, small_step) : vec3f_add(*hit_point, small_step);

    const vec3f_t ref_color = get_pixel_color(&ray_origin, &refracted_ray, spheres, lights, depth + 1);

    return vec3f_scale(ref_color, material->refraction_constant);
}

vec3f_t get_pixel_color(
    const vec3f_t *ray_origin, 
    const vec3f_t *ray_direction, 
    const array_t *spheres,
    const array_t *lights,
    const int depth)
{
    static const int MAX_RECURSION_DEPTH = 4;
    static const vec3f_t BACKGROUND_COLOR = { .4f, .5f, .8f };

    vec3f_t color_result = BACKGROUND_COLOR;
    material_t sphere_material;
    vec3f_t hit_point, normal;

    if (depth > MAX_RECURSION_DEPTH || 
        !cast_ray_to_spheres(ray_origin, ray_direction, spheres, &hit_point, &normal, &sphere_material))
        return color_result;

    float diffuse_intensity = calculate_diffuse_lighting(lights, &hit_point, &normal, &sphere_material, spheres);

    const vec3f_t diffuse_color = vec3f_scale(sphere_material.color, diffuse_intensity);
    const vec3f_t specular_color = 
        calculate_specular_lighting(lights, ray_direction, &hit_point, &normal, &sphere_material, spheres);
    const vec3f_t reflected_color = 
        calculate_reflected_color(ray_direction, &normal, &hit_point, spheres, lights, &sphere_material, depth);
    const vec3f_t refracted_color = 
        calculate_refracted_color(ray_direction, &normal, &hit_point, spheres, lights, &sphere_material, depth);
    
    color_result = vec3f_add(vec3f_add(diffuse_color, specular_color), vec3f_add(reflected_color, refracted_color));    

    // NOTE: Normalizing the color result makes cool cartoony effect
    return clamp_color(&color_result);
}

void draw_picture(vec3f_t *picture, const array_t *spheres, const array_t *lights)
{
    static const float fov = 90.0f; // Field of view
    static const float aspect_ratio = (float)WIDTH / (float)HEIGHT;
    
    const float screen_half_width = tanf(fov * M_PI / 360.0f);

    # pragma omp parallel for
    for (uint16_t y = 0; y < HEIGHT; y++)
        for (uint16_t x = 0; x < WIDTH; x++)
        {
            const float pixel_x = x + 0.5f - WIDTH / 2.0f;
            const float pixel_y = -(y + 0.5f - HEIGHT / 2.0f);
            const float pixel_z = -HEIGHT / (2 * tanf(fov * M_PI / 360.0f));
            const vec3f_t ray_direction = vec3f_normalize((vec3f_t){pixel_x, pixel_y, pixel_z});
            
            picture[y * WIDTH + x] = get_pixel_color(&(vec3f_t){0, 0, 0}, &ray_direction, spheres, lights, 0);
        }
}

void save_image(vec3f_t *picture)
{
    FILE *outfile = fopen("image.ppm", "w");    

    fprintf(outfile, "P6\n%d %d\n255\n", WIDTH, HEIGHT);
    for (uint32_t i = 0; i < WIDTH * HEIGHT; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            uint8_t pixel = 255 * fminf(((float *)(picture + i))[j], 1.f);

            fprintf(outfile, "%c", pixel);
        }
    }

    fclose(outfile);
}

void render(const array_t *spheres, const array_t *lights)
{
    vec3f_t *picture = malloc(WIDTH * HEIGHT * sizeof(vec3f_t));

    draw_picture(picture, spheres, lights);
    save_image(picture);
}

int main()
{
    material_t material_1 = {
        .color = {.7f, .3f, .2f},
        .shininess = 32.f,
        .refraction_index = 1.f, 
        .diffuse_constant = .8f,
        .specular_constant = .5f,
        .reflection_constant = .1,
        .refraction_constant = 0.f
    };

    material_t material_2 = {
        .color = {.9f, .1f, .2f},
        .shininess = 16.f,
        .refraction_index = 1.f,
        .diffuse_constant = 1.f,
        .specular_constant = 0.4f,
        .reflection_constant = 0.2f,
        .refraction_constant = 0.6f
    };

    material_t material_3 = {
        .color = {.3f, .3f, .5f},
        .shininess = 1.f,
        .refraction_index = 1.f,
        .diffuse_constant = .7f,
        .specular_constant = .1f,
        .reflection_constant = .9f,
        .refraction_constant = 0.f
    };

    material_t material_4 = {
        .color = {.3f, .7f, .1f},
        .shininess = 1.f,
        .refraction_index = 1.f,
        .diffuse_constant = .8f,
        .specular_constant = .0f,
        .reflection_constant = .1f,
        .refraction_constant = 0.f
    };

    sphere_t spheres[5];
    spheres[0] = (sphere_t){{0.f, 0.f, -16.f}, 2.f, material_1};
    spheres[1] = (sphere_t){{4.f, -4.f, -12.f}, 4.f, material_2};
    spheres[2] = (sphere_t){{6.f, 5.f, -9.f}, 2.5f, material_3};
    spheres[3] = (sphere_t){{4.f, -7.f, -24.f}, 1.f, material_1};
    spheres[4] = (sphere_t){{-40.f, 18.f, -54.f}, 9.f, material_4};
    array_t sphere_array = {spheres, sizeof(spheres) / sizeof(sphere_t)};

    light_t lights[3];
    lights[0] = (light_t){{-10.f, 10.f, -20.f}, 1.5f};
    lights[1] = (light_t){{30.f, 50.f, -25.f}, 1.2f};
    lights[2] = (light_t){{40.f, 0.f, 10.f}, 0.8f};
    array_t light_array = {lights, sizeof(lights) / sizeof(light_t)};

    render(&sphere_array, &light_array);
    return 0;
}
