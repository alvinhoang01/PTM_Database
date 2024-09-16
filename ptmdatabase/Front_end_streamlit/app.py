import streamlit as st
import psycopg2
from psycopg2 import sql
import bcrypt
from datetime import datetime
from streamlit_option_menu import option_menu

# --- Database connection setup ---
def connect_db():
    return psycopg2.connect(
        host="localhost",
        database="user_auth",
        user="postgres",
        password="Kaka2001!",
    )

# --- User sign-up function ---
def sign_up(email, password):
    try:
        # This will likely fail, but doing it for testing purposes
        salt = bcrypt.gensalt()
        hashed_password = bcrypt.hashpw(password, salt)  # Passing password as string (will fail)

        # Print hashed password details for debugging
        print(f"Hashed Password: {hashed_password}")

    except Exception as e:
        # Catch any errors during hashing and display them
        st.error(f"Error during password hashing: {e}")
        return

    # Database connection and recording
    try:
        conn = connect_db()
        cur = conn.cursor()

        # Insert the email and hashed password (stored in bytes)
        cur.execute(
            sql.SQL("INSERT INTO users (email, hashed_password, created_at) VALUES (%s, %s, %s)"),
            (email, hashed_password, datetime.now())
        )
        conn.commit()
        st.success(f"User {email} successfully registered!")

    except psycopg2.IntegrityError:
        conn.rollback()
        st.error(f"Email {email} already exists.")
    except Exception as db_error:
        st.error(f"Database error: {db_error}")
    finally:
        try:
            cur.close()
            conn.close()
        except Exception as cleanup_error:
            st.error(f"Error during connection close: {cleanup_error}")




# --- User login function ---
def login(email, password):
    conn = connect_db()
    cur = conn.cursor()

    cur.execute(
        sql.SQL("SELECT hashed_password FROM users WHERE email = %s"),
        (email,)
    )

    user_record = cur.fetchone()
    cur.close()
    conn.close()

    if user_record:
        stored_hashed_password = user_record[0]
        if bcrypt.checkpw(password.encode('utf-8'), stored_hashed_password.encode('utf-8')):
            return True
    return False

# --- Display Sidebar and Pages ---
def display_sidebar_and_pages(email):
    st.sidebar.title(f"Welcome, {email}!")
    st.sidebar.markdown("#### Select a page to navigate:")

    with st.sidebar:
        page = option_menu(
            menu_title=None,
            options=["Home Page", "Database Generation", "Matrix Analysis"],
            icons=["house", "database", "bar-chart-line"],
            menu_icon="cast",
            default_index=0,
            orientation="vertical",
            styles={
                "nav-link": {
                    "font-size": "16px",
                    "text-align": "left",
                    "margin": "5px",
                    "--hover-color": "#f0f0f0",
                },
                "icon": {
                    "font-size": "18px",
                    "color": "#008cba",
                },
                "container": {
                    "padding": "5px",
                    "background-color": "transparent"
                },
                "nav-link-selected": {
                    "background-color": "rgba(0, 123, 255, 0.15)",
                    "font-weight": "bold",
                    "color": "#008cba",
                },
            },
        )

    # --- Load Pages Based on Sidebar Selection ---
    if page == "Home Page":
        from Home_page import main as show_home_page
        show_home_page()

    elif page == "Database Generation":
        from Components.Database_Generation import main as show_database_page
        show_database_page()

    elif page == "Matrix Analysis":
        from Components.Matrix_analysis import main as show_matrix_analysis_page
        show_matrix_analysis_page()

    # Add logout button to the sidebar
    if st.sidebar.button('🔒 Logout'):
        st.session_state['authenticated'] = False
        st.session_state['email'] = ''
        st.session_state['signup'] = False
        # Since we are no longer using `experimental_rerun`, this should reload properly based on session state

# --- Streamlit interface for sign-up and login ---
def main():
    st.set_page_config(
        page_title="Proteoform Database Generation App",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="collapsed",
    )

    # --- Initialize session state ---
    if 'authenticated' not in st.session_state:
        st.session_state['authenticated'] = False
        st.session_state['email'] = ''
        st.session_state['signup'] = False

    if st.session_state['authenticated']:
        # If authenticated, show the main content
        display_sidebar_and_pages(st.session_state['email'])
    else:
        # If not authenticated, show login or sign-up form
        if st.session_state['signup']:
            # Sign-up form
            st.subheader("Create New Account")

            email = st.text_input("Email")
            password = st.text_input("Password", type='password')
            confirm_password = st.text_input("Confirm Password", type='password')

            sign_up_clicked = st.button("Sign Up")
            back_to_login_clicked = st.button("Back to Login")

            if sign_up_clicked:
                if password == confirm_password:
                    sign_up(email, password)
                    st.session_state['signup'] = False  # Return to login after sign-up
                    st.success("Please log in with your new credentials.")
                else:
                    st.error("Passwords do not match.")

            elif back_to_login_clicked:
                st.session_state['signup'] = False

        else:
            # Login form
            st.subheader("Login")

            email = st.text_input("Email")
            password = st.text_input("Password", type='password')

            login_clicked = st.button("Login")
            go_to_signup_clicked = st.button("Go to Sign Up")

            if login_clicked:
                if login(email, password):
                    st.success(f"Logged in successfully as {email}!")
                    st.session_state['authenticated'] = True
                    st.session_state['email'] = email
                else:
                    st.error("Email/password is incorrect")

            elif go_to_signup_clicked:
                st.session_state['signup'] = True

if __name__ == '__main__':
    main()
